# A collection of functions used to process nano-FTIR data and fitting
# Author: Yinming Shao (2026)

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from snippet import subplots, Sky, TOLO
from scipy.optimize import least_squares
from scipy.signal import savgol_filter, find_peaks
from scipy.special import hankel1
import copy

# ==========================================
# 1. 纯净版：只负责单文件读取的 readnanoFTIR
# ==========================================
def readnanoFTIR(file_path, ref=False, harmonics=(1,2,3,4,5)):
    """
    修改点：去掉了原本负责平均的 if isinstance(file_path, list): 逻辑。
    现在它是一个纯粹的数据读取器，返回原始的 info 和 dataframe。
    """
    if isinstance(file_path, list):
        raise ValueError("readnanoFTIR 不再支持传入列表！请使用 avg_scan 函数进行多文件平均。")
        
    info = {} 
    with open(file_path, 'r', encoding='utf-8') as file:
        next(file)
        for line in file:
            line = line.replace('0xa0', ' ').replace('\xa0', ' ')
            if line.startswith('#'):
                if ':' in line:
                    key, value = line[1:].split(':', 1)
                    info[key.strip()] = value.strip()

    df = pd.read_csv(file_path, sep='\t', comment='#')
    
    scansize = float(info['Scan Area (X, Y, Z)'].split('\t')[1]) 
    pixels = float(info['Pixel Area (X, Y, Z)'].split('\t')[1])
    pixelsize = scansize*1e3/pixels
    info['pixelsize_nm'] = round(pixelsize,3)

    if ref:
        ref_info = {}
        with open(ref, 'r', encoding='utf-8') as file:
            next(file)
            for line in file:
                line = line.replace('0xa0', ' ')
                if line.startswith('#') and ':' in line:
                    key, value = line[1:].split(':', 1)
                    ref_info[key.strip()] = value.strip()

        ref_df = pd.read_csv(ref, sep='\t', comment='#')
        df = pd.merge(df, ref_df, on='Wavenumber', how='left', suffixes=('', '_ref'))
        
        for n in harmonics:
            A = f'O{n}A'
            P = f'O{n}P'
            Aref = f'{A}_ref'
            Pref = f'{P}_ref'

            if A in df.columns and P in df.columns and Aref in df.columns and Pref in df.columns:
                df.loc[df[Aref] == 0, Aref] = np.nan
                df[A] = df[A] / df[Aref]
                df[P] = df[P] - df[Pref]
                df[P] = (df[P] + np.pi) % (2*np.pi) - np.pi

        ref_cols = [c for c in df.columns if c.endswith('_ref')]
        df = df.drop(columns=ref_cols)
        info['external_reference_applied'] = True
        info['external_reference_file'] = ref
        if 'Description' in ref_info:
            info['external_reference_description'] = ref_info['Description']
    else:
        info['external_reference_applied'] = False

    harmonics_set = set(harmonics)
    drop_cols = []
    for c in df.columns:
        if len(c) >= 3 and c[0] == 'O' and c[-1] in ['A', 'P'] and c[1:-1].isdigit():
            n = int(c[1:-1])
            if n not in harmonics_set:
                drop_cols.append(c)
    if len(drop_cols) > 0:
        df = df.drop(columns=drop_cols)

    for k, v in info.items():
        print(f"{k}: {v}")

        
    return {'info': info, 'data': df}

# ==========================================
# 2. 辅助运算函数 (供 avg_scan 内部调用)
# ==========================================
def _normalize_raw_df(df_data, normby, harmonics):
    """用于实现 '先归一化再平均' 的长表矩阵处理"""
    df_norm = df_data.copy()
    for n in harmonics:
        A = f'O{n}A'
        P = f'O{n}P'
        if A not in df_norm.columns or P not in df_norm.columns:
            continue
            
        # 提取基底参考平均
        if len(normby) == 2:
            mask = (df_norm['Column'] >= normby[0]) & (df_norm['Column'] < normby[1])
            ref = df_norm[mask].groupby('Wavenumber')[[A, P]].mean()
        elif len(normby) == 4:
            mask1 = (df_norm['Column'] >= normby[0]) & (df_norm['Column'] < normby[1])
            mask2 = (df_norm['Column'] >= normby[2]) & (df_norm['Column'] < normby[3])
            ref1 = df_norm[mask1].groupby('Wavenumber')[[A, P]].mean()
            ref2 = df_norm[mask2].groupby('Wavenumber')[[A, P]].mean()
            ref = (ref1 + ref2) / 2
        else:
            continue
            
        # 广播相除/相减
        df_norm = df_norm.merge(ref, on='Wavenumber', suffixes=('', '_ref'))
        df_norm[A] = df_norm[A] / df_norm[f'{A}_ref']
        df_norm[P] = df_norm[P] - df_norm[f'{P}_ref']
        df_norm[P] = (df_norm[P] + np.pi) % (2*np.pi) - np.pi
        df_norm.drop(columns=[f'{A}_ref', f'{P}_ref'], inplace=True)
    return df_norm

def _to_complex(df, harmonics):
    """极坐标(A,P)转复平面(Re,Im)"""
    df_c = df.copy()
    for n in harmonics:
        A, P = f'O{n}A', f'O{n}P'
        if A in df_c.columns and P in df_c.columns:
            df_c[f'O{n}_Re'] = df_c[A] * np.cos(df_c[P])
            df_c[f'O{n}_Im'] = df_c[A] * np.sin(df_c[P])
            df_c.drop(columns=[A, P], inplace=True)
    return df_c

def _to_polar(df, harmonics):
    """复平面(Re,Im)转回极坐标(A,P)"""
    df_p = df.copy()
    for n in harmonics:
        Re, Im = f'O{n}_Re', f'O{n}_Im'
        if Re in df_p.columns and Im in df_p.columns:
            df_p[f'O{n}A'] = np.sqrt(df_p[Re]**2 + df_p[Im]**2)
            df_p[f'O{n}P'] = np.arctan2(df_p[Im], df_p[Re])
            df_p.drop(columns=[Re, Im], inplace=True)
    return df_p

# ==========================================
# 3. 核心大管家：独立出来的 avg_scan
# ==========================================
def avg_scan(registry, label_main, label_repeat,
             normby_main, normby_repeat,
             ch_params,
             wran=(195, 400), avg_rows=2, pointspec=(5, 75),
             harmonics=(1, 2, 3, 4, 5), ref=False, 
             avg_mode='complex', norm_first=False):
    """
    全新升级版：整合所有对齐、裁剪、标量/复数平均逻辑。
    新增参数 norm_first=True 可以控制先执行背景归一化，再进行平均。
    """
    print(f"\n=============================================")
    print(f"🔧 [avg_scan] 开始处理: {label_main}")
    print(f"   平均模式: {avg_mode.upper()} | 运算顺序: {'先归一化，再平均' if norm_first else '先平均，再归一化'}")
    
    # 1. 独立读取两份数据
    res_main = readnanoFTIR(registry[label_main]["NearField"]["path"], ref=ref, harmonics=harmonics)
    res_repeat = readnanoFTIR(registry[label_repeat]["NearField"]["path"], ref=ref, harmonics=harmonics)
    
    df_main = res_main['data'].copy()
    df_repeat = res_repeat['data'].copy()

    # 2. 如果开启先归一化，此时独立处理各自背景
    if norm_first:
        df_main = _normalize_raw_df(df_main, normby_main, harmonics)
        df_repeat = _normalize_raw_df(df_repeat, normby_repeat, harmonics)

    # 3. 计算 Drift 和平移量
    drift_idx = 1 if len(normby_main) == 2 else 2
    drift = normby_repeat[drift_idx] - normby_main[drift_idx]
    column_shifts = [0, -drift]
    
    df_main['Column'] = df_main['Column'] + column_shifts[0]
    df_repeat['Column'] = df_repeat['Column'] + column_shifts[1]

    # 4. 求交集并对齐裁剪
    common_cols = sorted(set(df_main['Column'].unique()) & set(df_repeat['Column'].unique()))
    df_main = df_main[df_main['Column'].isin(common_cols)]
    df_repeat = df_repeat[df_repeat['Column'].isin(common_cols)]

    # 5. 执行核心平均逻辑 (解耦完成)
    if avg_mode == 'complex':
        df_main_c = _to_complex(df_main, harmonics)
        df_repeat_c = _to_complex(df_repeat, harmonics)
        df_mean = pd.concat([df_main_c, df_repeat_c]).groupby(['Wavenumber', 'Column'], as_index=False).mean()
        df_mean = _to_polar(df_mean, harmonics)
    elif avg_mode == 'scalar':
        df_mean = pd.concat([df_main, df_repeat]).groupby(['Wavenumber', 'Column'], as_index=False).mean()
    else:
        raise ValueError(f"avg_mode 必须是 'complex' 或 'scalar', 收到的是 '{avg_mode}'")

    # 6. 重置列索引 (确保从 0 开始)
    col_map = {c: i for i, c in enumerate(common_cols)}
    df_mean['Column'] = df_mean['Column'].map(col_map)
    n_px = df_mean['Column'].max() + 1
    
    # 7. 更新 Info
    info = res_main['info'].copy()
    orig_px = int(info['Pixel Area (X, Y, Z)'].split('\t')[1])
    orig_scan = float(info['Scan Area (X, Y, Z)'].split('\t')[1])
    px_size_um = orig_scan / orig_px

    pa_parts = info['Pixel Area (X, Y, Z)'].split('\t')
    pa_parts[1] = str(n_px)
    info['Pixel Area (X, Y, Z)'] = '\t'.join(pa_parts)
    
    sa_parts = info['Scan Area (X, Y, Z)'].split('\t')
    sa_parts[1] = f'{px_size_um * n_px:.4f}'
    info['Scan Area (X, Y, Z)'] = '\t'.join(sa_parts)
    info['pixelsize_nm'] = round(px_size_um * 1e3, 3)
    
    # 8. 计算新的 normby 区间并保护
    left_crop = max(0, -drift)
    normby_avg = [max(0, min(x - left_crop, n_px)) for x in normby_main]
    
    # 🛑 核心修复：防止右侧被裁空导致 NaN 崩溃
    if len(normby_avg) == 4 and normby_avg[2] >= normby_avg[3]:
        print(f"   ⚠️ 警告: 漂移导致右侧基底完全消失！为防止报错，自动退化为左侧单基底: {normby_avg[:2]}")
        normby_avg = normby_avg[:2]
        
    # 打印你关心的计算过程
    print(f"   🔍 漂移计算: drift = {drift} px  -> column_shifts = {column_shifts}")
    print(f"   🔍 裁剪统计: 原始 {orig_px} 列, 重合保留 {n_px} 列 (左侧裁掉 {left_crop} 列)")
    print(f"   ✅ [关键输出] 修正后的基底区间 (normby_avg): {normby_main} -> {normby_avg}")
    print(f"=============================================\n")

    # 9. 封装最终数据
    df = {'info': info, 'data': df_mean}
    nf_cfg = copy.deepcopy(registry[label_main]["NearField"])
    nf_cfg["wran"] = list(wran)
    nf_cfg["avg"] = avg_rows
    nf_cfg["pointspec"] = list(pointspec)
    nf_cfg["normby"] = normby_avg
    nf_cfg["channels"] = copy.deepcopy(ch_params)

    # 10. 传给作图函数
    for ch, n in zip(["O1", "O2"], [1, 2]):
        if ch not in ch_params: continue
        ch_cfg = ch_params[ch]
        
        # 如果已经先归一化过了，plotnanoFTIR 这里传相同的区间，除以 1、减去 0，等价于跳过二次归一化，但保留了画线功能！
        fig, imgdata = plotnanoFTIR(
            df, n=n,
            wran=list(wran),
            normby=normby_avg,
            amp_range=ch_cfg.get("v_range_amp", [0, 2]),
            phase_range=ch_cfg.get("v_range_phase", [-1, 1]),
            linecut=ch_cfg.get("linecut_freq", False),
            avg=avg_rows,
            pointspec=list(pointspec),
        )
        df[ch] = imgdata

    return df, nf_cfg

# 【注意】：下面的 plotnanoFTIR 函数保持你原有的即可，无需任何修改！

# 从这里往下保留你原本的 plotnanoFTIR 等函数...

def plotnanoFTIR(df,n=1,wran=[180,1000],figAx=False,normby=False,amp_range=[0,2],
                 phase_range=[-np.pi,np.pi],shading='auto',colorbar=True,linecut=False,
                 avg=0,subtract_bg=False,pointspec=True,savefig=False,rasterized=False):
    """
    plot nano-FTIR hyperspectra linescan from Neaspec and return O1, O2 dataframe
    
    n: 1 (Default plot the 1st harmonic amplitude and phase: O1A and O1P)
    Default colormap is Sky for amplitude and RdBu for phase
    wran: frequency range in wavenumebrs, default [180,1000] 
    
    normby: if not False, normalize each spectra by some references
    'avg': normalize by the average of all spectra, good for enhancing dispersive features
    'substrate': normalize by the substrate lines, by default use the first 5 lines as the
    substrate. Alternatively, specify the line number range as substrate: 
    e.g. [0,20], where line number 0 is the start of the hyperspectra line scan (usually on substrate)
    or [0,10,70,80], where both the start [0,10] and end [70,80] are substrate
    
    linecut: if not False, will add a linecut and export the data at the specified wavenumber.
    e.g. linecut=330 will add a dashed line (red) at 330 cm-1 and include the linecut into the
    imgdata dict ({'wavenumber':ff,'distance_um':ll,'O'+str(n)+'A':amp,'O'+str(n)+'P':phase})
    For linecuts at more than one freq, use linecut = [freq1, freq2,freq3,...]
    
    avg: if not 0, will average the specified number of rows around the freq linecut to improve
    the signal-to-noise of the linecut. The highest frequneyc resolution of Neaspec is 3.3cm-1,
    limited by the range of interferometer (L=1.5 mm): 1/2L = 1/3mm = 3.3 cm-1.
    But the data from Neaspec is oversampled and have a nominal resolution of 0.8333 cm-1, 
    so an average of 4 freqs near the designated freqcut is a good starting point.
    
    subtract_bg: if not False, will subtract a linear background based on fitting at each freq
    If normby is True, add 1 to the subtracted amp data since the orginal normed value is 1.
    
    If both linecut and avg are not False, will plot a side panel on the left to display the
    averaged linecut as a function of distance.
    
    pointspec: default True, will average the entire sample region to extract non-dipsersive feature's
    point spectra. The sample region will be detemined by normby.
    if provided a specific range, will only average within the range specifled, e.g. [10,20]
    
    returns fig, imgdata
    
    """
    # process multiple linecuts
    def _parse_linecuts(linecut):
        """Return [] / [scalar] / list of scalars."""
        if linecut is False or linecut is None:
            return []
        if np.isscalar(linecut):
            return [float(linecut)]
        return [float(x) for x in linecut]
    
    
    if figAx:
        fig,ax1,ax2 = figAx
    else:
        fig,(ax1,ax2) = subplots(2,size=(8,8))
        
    px = max(df['data']['Column'])+1  # number of pixels
    pxinfo = int(df['info']['Pixel Area (X, Y, Z)'].split('\t')[1])
    scanlength = float(df['info']['Scan Area (X, Y, Z)'].split('\t')[1])
    if px != pxinfo:
        print('Inconsistent number of pixel in data and info: %d, %d!'%(px,pxinfo))
    else:
        amp = df['data'].pivot(index='Wavenumber',columns='Column',values='O%dA'%n)
        phase = df['data'].pivot(index='Wavenumber',columns='Column',values='O%dP'%n)
        # only keep the data that are within the wran to avoid lines with NaNs
        amp = amp.query('@wran[0]<Wavenumber<@wran[1]')
        phase = phase.query('@wran[0]<Wavenumber<@wran[1]')
        freq = amp.index
        length = amp.columns*scanlength/px
        ll,ff = np.meshgrid(length, freq)

        # make sure normby range is valid and exceeding total scan size
        edge = None   # default: no normby edge markers
        if normby:
            if normby == 'avg':
                amp_avg = amp.mean(axis=1) # avreage of all columns/distances
                phase_avg = phase.mean(axis=1)
            elif type(normby) == list:
                # only use lines specfied by normby as substrate (reference)
                if any(x > px for x in normby):
                    print('normby contains indices larger than the total scan size px = %d!'%px)
                if len(normby) == 2:
                    # [start,end] of the substrate on one side of the image
                    amp_avg = amp.iloc[:,normby[0]:normby[1]].mean(axis=1)
                    phase_avg = phase.iloc[:,normby[0]:normby[1]].mean(axis=1)
                    if normby[0]>normby[1]:
                        # for cases where last lines are substrate, e.g. [80,-1]
                        edge = normby[0]*df['info']['pixelsize_nm']/1e3
                    else: # normal case where first lines are substrate, e.g. [0,10]
                        edge = normby[1]*df['info']['pixelsize_nm']/1e3
                    ax1.axhline(edge,ls='--',c='silver')
                    ax2.axhline(edge,ls='--',c='silver')
                elif len(normby) == 4:
                    # [start1,end1,start2,end2] for substrate on both side of the image
                    amp_avg_1 = amp.iloc[:, normby[0]:normby[1]].mean(axis=1)
                    amp_avg_2 = amp.iloc[:, normby[2]:normby[3]].mean(axis=1)
                    amp_avg = (amp_avg_1 + amp_avg_2) / 2

                    phase_avg_1 = phase.iloc[:, normby[0]:normby[1]].mean(axis=1)
                    phase_avg_2 = phase.iloc[:, normby[2]:normby[3]].mean(axis=1)
                    phase_avg = (phase_avg_1 + phase_avg_2) / 2
                    edge1 = normby[1]*df['info']['pixelsize_nm']/1e3
                    edge2 = normby[2]*df['info']['pixelsize_nm']/1e3
                    for edge in [edge1, edge2]:
                        ax1.axhline(edge, ls='--', c='silver')
                        ax2.axhline(edge, ls='--', c='silver')
                    edge = (edge1, edge2)  # tuple of edges
            
            # more robust way to do the normalizations
            amp = amp.div(amp_avg, axis=0)
            phase = phase.sub(phase_avg, axis=0)
            #for i in amp.columns:
                #amp.iloc[:,i] = amp.iloc[:,i]/amp_avg
                #phase.iloc[:,i] = phase.iloc[:,i] - phase_avg
               
        if subtract_bg:
            # For each frequency "row" in amp, fit amplitude vs. distance to a line
            for i in range(len(amp)):
                p = np.polyfit(length, amp.iloc[i], 1)
                amp.iloc[i] -= np.polyval(p, length) - 1
            # Do the same for phase
            for i in range(len(phase)):
                p = np.polyfit(length, phase.iloc[i], 1)
                phase.iloc[i] -= np.polyval(p, length)
            
        # 2D colorplot
        for ax,data,cmap,title,vran in zip([ax1,ax2],[amp,phase],[Sky,'RdBu'],
                                      [r'$\rm{S_%d}$ (norm.)'%(n),r'$\rm{\phi_%d}$ (norm.)'%(n)],[amp_range,phase_range]):
            img = ax.pcolormesh(ff, ll, data,cmap=cmap,shading=shading,vmin=vran[0],vmax=vran[1],rasterized=rasterized)
            axt = ax.twiny()
            axt.set(xlabel='Energy (meV)',xlim=np.array(wran)/8.065)
            if colorbar:
                cb = fig.colorbar(img,ax=ax,shrink=0.8)
                cb.ax.set_title(title,rotation=90,x=-1.5,y=0.35)
            ax.set(ylabel=r'Distance ($\rm{\mu m}$)',xlabel=r'Frequency (cm$^{-1}$)',xlim=wran)

        # allow mutiple linecuts to be specified now
        if linecut:
            linecuts = _parse_linecuts(linecut)
            lp = pd.DataFrame({'distance_um': np.asarray(length)})

            for lc in linecuts:
                freq_index = amp.index.get_indexer([lc], method='nearest')[0]
                f_actual = float(amp.index[freq_index])
                ax1.axvline(lc, ls='--', c='r', alpha=0.5)
                ax2.axvline(lc, ls='--', c='r', alpha=0.5)
                if avg:
                    # average a few freqs around the freq_index defined by the number of avg, axis=0 is crucial
                    i0 = max(0, freq_index - round(avg/2))
                    i1 = min(len(amp), freq_index + round(avg/2) + 1)
                    amp_lc = amp.iloc[i0:i1].mean(axis=0).to_numpy()
                    phase_lc = phase.iloc[i0:i1].mean(axis=0).to_numpy()
                else:
                    amp_lc = amp.iloc[freq_index].to_numpy()
                    phase_lc = phase.iloc[freq_index].to_numpy()

                # use actual frequency in column label so you know what was really used
                tag = f'{"%.2f" % f_actual}cm-1'
                lp[f'{tag}_O{n}A'] = amp_lc
                lp[f'{tag}_O{n}P'] = phase_lc
                

        if linecut and avg:
            # allow scalar or list
            lcs = [linecut] if np.isscalar(linecut) else list(linecut)

            if len(lcs) == 1:
                # original style: one panel with amp + phase
                lc = lcs[0]
                freq_index = amp.index.get_indexer([lc], method='nearest')[0]
                i0 = max(0, freq_index - round(avg/2))
                i1 = min(len(amp), freq_index + round(avg/2) + 1)

                figlc, ax = subplots(size=(7,3))
                # averaged curves (stored in lp)
                ax.plot(lp['distance_um'], lp.iloc[:,1], lp['distance_um'], lp.iloc[:,2])
                ax.legend(['Amplitude','Phase'], title=r'$\rm{\omega=%g\,cm^{-1}}$'%lc)
                
                # gray curves used in averaging
                for i in range(i0, i1):
                    ax.plot(lp['distance_um'], amp.iloc[i], c='gray', alpha=0.2, zorder=0)
                    ax.plot(lp['distance_um'], phase.iloc[i], c='gray', alpha=0.2, zorder=0)

                if edge:
                    if isinstance(edge, (list, tuple)):
                        ax.axvspan(0, edge[0], color='gray', alpha=0.1)
                        ax.axvspan(edge[1], np.array(lp['distance_um'])[-1], color='gray', alpha=0.1)
                    else:
                        ax.axvspan(0, edge, color='gray', alpha=0.1)

                ax.set(xlabel=r'Distance ($\rm{\mu m}$)',
                       ylabel=r'Averaged $S_%d,\,\phi_%d$'%(n,n), xlim=0)

            else:
                # multiple linecuts: amplitude and phase in separate panels
                figlc, (axA, axP) = subplots(1, 2, size=(10,5))

                for j, lc in enumerate(lcs):
                    freq_index = amp.index.get_indexer([lc], method='nearest')[0]
                    i0 = max(0, freq_index - round(avg/2))
                    i1 = min(len(amp), freq_index + round(avg/2) + 1)

                    # gray curves used in averaging window
                    for i in range(i0, i1):
                        axA.plot(lp['distance_um'], amp.iloc[i], c='gray', alpha=0.2, zorder=0)
                        axP.plot(lp['distance_um'], phase.iloc[i], c='gray', alpha=0.2, zorder=0)

                    # averaged curves from lp (assumes columns are [dist, lc1A, lc1P, lc2A, lc2P, ...])
                    axA.plot(lp['distance_um'], lp.iloc[:, 1 + 2*j], label=f'{lc:g}')
                    axP.plot(lp['distance_um'], lp.iloc[:, 2 + 2*j], label=f'{lc:g}')

                axA.legend(title=r'$\omega$ (cm$^{-1}$)', ncol=2, fontsize=9)

                if edge:
                    for ax_ in (axA, axP):
                        if isinstance(edge, (list, tuple)):
                            ax_.axvspan(0, edge[0], color='gray', alpha=0.1)
                            ax_.axvspan(edge[1], np.array(lp['distance_um'])[-1], color='gray', alpha=0.1)
                        else:
                            ax_.axvspan(0, edge, color='gray', alpha=0.1)

                axA.set(xlabel=r'Distance ($\rm{\mu m}$)', ylabel=rf'Averaged $S_{n}$', xlim=0)
                axP.set(xlabel=r'Distance ($\rm{\mu m}$)', ylabel=rf'Averaged $\phi_{n}$', xlim=0)
            
        if pointspec:
            if type(pointspec) == list:
                # this will also add some line indicators on the colorplot to indicate the regions used for average
                ampspec = amp.iloc[:, pointspec[0]:pointspec[1]].mean(axis=1) 
                phaspec = phase.iloc[:, pointspec[0]:pointspec[1]].mean(axis=1)
            else:
                # assume the first lines are subtrate, and lines after normby[1] are all sample
                ampspec = amp.iloc[:, normby[1]:].mean(axis=1) # average all sample region and generate a single spectra
                phaspec = phase.iloc[:, normby[1]:].mean(axis=1)
            spec = pd.concat([ampspec,phaspec],axis=1,keys=['O'+str(n)+'A','O'+str(n)+'P']).reset_index()
                
        imgdata = {'wavenumber':ff,'distance_um':ll,'linecut':lp,'pointspectra':spec,
                   'O'+str(n)+'A':amp,'O'+str(n)+'P':phase, 'normby_px': normby}     
       
        return fig,imgdata   
  

## Example Usage ##  
#df = readnanoFTIR('nanoFTIR/NSLS2/Day1/300nm/2024-04-11 185547 NF LS CSB_300nm_baxis.txt')
# f,imgdata = plotnanoFTIR(df,n=1,wran=[195,400],normby=[0,29],
#                          amp_range=[0,2.2],phase_range=[-2,2],linecut=[335,338],avg=3)
# df # imgdata['linecut']




def fit_cavity_prefactor_compare(y_um, s, xr=(0.1, 3.7), yc_um=1.9,
                                 prefactor='none',   # 'none' | '1/sqrtx' | '1/x' | 'powerlaw' | 'hankel'
                                 R_nm=25.0,
                                 a0=1.0, fit_a=False,   # used only for 'powerlaw'
                                 fit_yc=False, yc_bounds=None,
                                 lam_bounds=None, lam_bound_factor=(0.7, 1.2),
                                 win=7, prom=0.007, robust=True, lam0_guess=None, edges='double'):
    """
    Minimal cavity fit with symmetric left/right edge distances.

    Non-Hankel models:
      s(y) = B + A * sin(2q(y-yc)+phi) * cosh(alpha_env*(y-yc)) * G(y)

      Here alpha_env is an *effective envelope* parameter (units 1/um).
      It is not necessarily the intrinsic Im(q_p); it absorbs multiple effects.

    Hankel model (full cylindrical waves from both edges, symmetric):
      s(y) = B + A * Re{ exp(i phi) [H0^(1)(2 q_p xL) + H0^(1)(2 q_p xR)] }

      where q_p = q + i*q_imag, q = 2π/lambda.
      IMPORTANT: we use 2*q_p*x (round-trip distance) to represent a tip-launched
      tip-edge-tip pathway. In the asymptotic form this gives exp(-2*q_imag*x)
      attenuation versus distance-from-edge x.

    Units:
      y, lambda in um
      q in rad/um internally, q_cm^-1 returned
      alpha_env or q_imag in 1/um (and cm^-1 returned)
    """
    import numpy as np
    from scipy.optimize import least_squares
    from scipy.signal import savgol_filter, find_peaks
    from scipy.special import hankel1

    # -------------------------
    # Crop to fit window
    # -------------------------
    m = (y_um >= xr[0]) & (y_um <= xr[1])
    y = np.asarray(y_um[m], float)
    ss = np.asarray(s[m], float)

    # -------------------------
    # Distances to the two cavity edges (symmetric construction)
    # Note: here edges are tied to xr. Later you may want separate edges_um.
    # -------------------------
    R_um = R_nm * 1e-3
    eps = max(R_um, 1e-6)
    xL = np.clip(y - xr[0], eps, None)
    xR = np.clip(xr[1] - y, eps, None)

    # -------------------------
    # Initial lambda guess from adjacent maxima spacing
    # (standing-wave adjacent maxima spacing ≈ lambda/2)
    # Use a small smoothing window to avoid merging peaks.
    # -------------------------
    wwin_guess = 5  # odd
    prom_guess = 0.005
    ys_guess = savgol_filter(ss, wwin_guess, 2)
    pks, _ = find_peaks(ys_guess, prominence=prom_guess)

    if lam0_guess is not None:
        lam0 = float(lam0_guess)
    else:
        dmax0 = np.median(np.diff(y[pks])) if len(pks) >= 2 else 0.4
        lam0 = float(np.clip(2*dmax0, 0.05, 20.0))

    # Automatic lambda bounds around lam0 unless user provides lam_bounds
    if lam_bounds is None:
        lo = max(lam_bound_factor[0] * lam0, 1e-4)
        hi = max(lam_bound_factor[1] * lam0, lo + 1e-6)
        lam_bounds_use = (lo, hi)
    else:
        lam_bounds_use = lam_bounds

    # -------------------------
    # Initial guesses
    # -------------------------
    B0 = np.nanmedian(ss)
    A0 = 0.5 * (np.nanpercentile(ss, 95) - np.nanpercentile(ss, 5))
    alpha0 = 0.5  # used as alpha_env for non-Hankel; used as q_imag for Hankel
    phi0 = 0.0

    # yc fitting is only used in the non-Hankel branch (yc not present in current Hankel model)
    use_yc = bool(fit_yc) and (prefactor != 'hankel')
    if yc_bounds is None:
        yc_bounds = (yc_um - 0.2, yc_um + 0.2)

    # -------------------------
    # Parameter vector
    # Non-Hankel: alpha0 = alpha_env
    # Hankel:     alpha0 = q_imag  (Im(q_p))
    # -------------------------
    use_a = (prefactor == 'powerlaw')
    p0 = [B0, A0, alpha0, np.clip(lam0, lam_bounds_use[0], lam_bounds_use[1]), phi0]
    lb = [-np.inf, -np.inf, 0.0, lam_bounds_use[0], -np.pi]
    ub = [ np.inf,  np.inf, 50.0, lam_bounds_use[1],  np.pi]

    if use_yc:
        p0.append(float(yc_um))
        lb.append(float(yc_bounds[0]))
        ub.append(float(yc_bounds[1]))

    if use_a:
        p0.append(float(a0))
        if fit_a:
            lb.append(0.0); ub.append(4.0)
        else:
            lb.append(float(a0)); ub.append(float(a0))

    p0, lb, ub = map(np.array, (p0, lb, ub))

    def geom_pref(alpha_or_qimag, lam, a=None):
        """
        Returns the symmetric geometry factor.

        For prefactor='hankel', this returns a complex field:
          H0^(1)(2*q_p*xL) + H0^(1)(2*q_p*xR)
        where q_p = q + i*q_imag and q_imag = alpha_or_qimag.

        For other prefactors, returns a real-valued geometric weighting G(y).
        """
        if prefactor == 'none':
            g = np.ones_like(y)

        elif prefactor == '1/sqrtx':
            g = 1/np.sqrt(xL + R_um)
            if edges == 'double': g += 1/np.sqrt(xR + R_um)

        elif prefactor == '1/x':
            g = 1/(xL + R_um)
            if edges == 'double': g += 1/(xR + R_um)

        elif prefactor == 'powerlaw':
            g = 1/(xL**a + R_um**a)
            if edges == 'double': g += 1/(xR**a + R_um**a)

        elif prefactor == 'hankel':
            q = 2*np.pi / lam
            q_imag = alpha_or_qimag
            q_p = q + 1j*q_imag

            # round-trip distance: 2*q_p*x
            zL = 2 * q_p * xL
            zR = 2 * q_p * xR
            field = hankel1(0, zL)
            if edges == 'double': field += hankel1(0, zR)
            return field

        else:
            raise ValueError("prefactor must be 'none', '1/sqrtx', '1/x', 'powerlaw', or 'hankel'")

        g = g / max(np.nanmax(np.abs(g)), 1e-12)
        return g

    def unpack_theta(theta):
        B, A, alpha_or_qimag, lam, phi = theta[:5]
        idx = 5

        yc_fit = yc_um
        if use_yc:
            yc_fit = theta[idx]
            idx += 1

        a = None
        if use_a:
            a = theta[idx]

        return B, A, alpha_or_qimag, lam, phi, yc_fit, a

    def model(theta):
        B, A, alpha_or_qimag, lam, phi, yc_fit, a = unpack_theta(theta)

        if prefactor == 'hankel':
            # alpha_or_qimag is q_imag here (Im(q_p))
            hc = geom_pref(alpha_or_qimag, lam, a=None)
            h = np.real(np.exp(1j*phi) * hc)
            h = h / max(np.nanmax(np.abs(h)), 1e-12)
            return B + A * h

        # Non-Hankel: alpha_or_qimag is alpha_env (effective envelope parameter)
        alpha_env = alpha_or_qimag
        q = 2*np.pi / lam
        Fint = np.cosh(alpha_env * (y - yc_fit))
        Fint = Fint / max(np.nanmax(Fint), 1e-12)

        G = geom_pref(alpha_env, lam, a)
        return B + A * np.sin(2*q*(y - yc_fit) + phi) * Fint * G

    def resid(theta):
        return model(theta) - ss

    res = least_squares(
        resid, p0, bounds=(lb, ub),
        loss='soft_l1' if robust else 'linear',
        f_scale=max(np.std(ss)*0.3, 1e-6),
        max_nfev=30000
    )

    fit = model(res.x)
    B, A, alpha_or_qimag, lam_um, phi, yc_fit, a_fit = unpack_theta(res.x)

    q_um = 2*np.pi / lam_um
    n = len(ss)
    k = len(res.x)
    rss = np.sum((fit - ss)**2)
    rmse = np.sqrt(rss / max(n-k, 1))
    aic = n*np.log(max(rss/n, 1e-30)) + 2*k

    # Put the damping parameter into a consistent field depending on the branch
    params = {
        'B': B, 'A': A, 'lambda_p_um': lam_um, 'phi': phi,
        'a': a_fit, 'yc_um': yc_fit, 'prefactor': prefactor, 'R_nm': R_nm,
        'fit_yc': use_yc, 'lam0_um': lam0, 'lam_bounds_um': lam_bounds_use
    }
    derived = {
        'q_rad_per_um': q_um,
        'q_cm^-1': q_um * 1e4,
        'adjacent_max_spacing_um': lam_um / 2
    }

    if prefactor == 'hankel':
        q_imag_um = alpha_or_qimag
        params['q_imag_um^-1'] = q_imag_um
        derived['q_imag_cm^-1'] = q_imag_um * 1e4
        # amplitude decay vs distance-from-edge x behaves like exp(-2*q_imag*x) due to round trip
        derived['amp_decay_length_um'] = 1.0 / max(2*q_imag_um, 1e-12)
    else:
        alpha_env_um = alpha_or_qimag
        params['alpha_env_um^-1'] = alpha_env_um
        derived['alpha_env_cm^-1'] = alpha_env_um * 1e4

    # Dense grid purely for a smooth plotted curve -- the fit itself (RMSE/AIC/params
    # above) is unchanged, still computed on the actual data grid `y`. Real nano-FTIR
    # line profiles are often sparse enough relative to a short lambda_p that connecting
    # only the data points with straight segments makes the curve look polygonal/jagged.
    x_dense = np.linspace(xr[0], xr[1], 600)
    xL_d = np.clip(x_dense - xr[0], eps, None)
    xR_d = np.clip(xr[1] - x_dense, eps, None)

    if prefactor == 'hankel':
        q_p_d = (2*np.pi/lam_um) + 1j*alpha_or_qimag
        field_d = hankel1(0, 2*q_p_d*xL_d)
        if edges == 'double':
            field_d = field_d + hankel1(0, 2*q_p_d*xR_d)
        h_d = np.real(np.exp(1j*phi) * field_d)
        h_d = h_d / max(np.nanmax(np.abs(h_d)), 1e-12)
        y_dense = B + A * h_d
    else:
        alpha_env = alpha_or_qimag
        q_d = 2*np.pi / lam_um
        Fint_d = np.cosh(alpha_env * (x_dense - yc_fit))
        Fint_d = Fint_d / max(np.nanmax(Fint_d), 1e-12)

        if prefactor == 'none':
            g_d = np.ones_like(x_dense)
        elif prefactor == '1/sqrtx':
            g_d = 1/np.sqrt(xL_d + R_um)
            if edges == 'double': g_d = g_d + 1/np.sqrt(xR_d + R_um)
        elif prefactor == '1/x':
            g_d = 1/(xL_d + R_um)
            if edges == 'double': g_d = g_d + 1/(xR_d + R_um)
        else:  # 'powerlaw'
            g_d = 1/(xL_d**a_fit + R_um**a_fit)
            if edges == 'double': g_d = g_d + 1/(xR_d**a_fit + R_um**a_fit)
        g_d = g_d / max(np.nanmax(np.abs(g_d)), 1e-12)
        y_dense = B + A * np.sin(2*q_d*(x_dense - yc_fit) + phi) * Fint_d * g_d

    return {
        'mask': m,
        'x_fit_um': y,
        'y_fit': fit,
        'x_dense_um': x_dense,
        'y_dense': y_dense,
        'params': params,
        'derived': derived,
        'metrics': {'rss': rss, 'rmse': rmse, 'aic': aic},
        'success': res.success,
        'fit_object': res
    }


def compare_cavity_models(amplp, col,
                          xr=(0.2, 3.6), yc_um=1.9, fit_yc=True,
                          R_nm=25.0, a0=1.0,
                          win=7, prom=0.02,
                          prefactors=('none', '1/sqrtx', '1/x', 'powerlaw', 'hankel'),
                          ylim=(0.5, 0.9), xlim=None,
                          figsize=(6.2, 5.2), show_text=True, lam0_guess=None, edges='double'):
    """
    Fit and compare multiple cavity models for one line-profile column.

    Parameters
    ----------
    amplp : DataFrame
        Must contain 'distance_um' and the selected column.
    col : str
        Example: '340.83cm-1_O2A'
    """
    x = amplp['distance_um'].values
    y = amplp[col].values

    # Default the displayed x-range to the fit window itself, rather than the full
    # data range -- otherwise the oscillation the fit actually targets gets squeezed
    # into a sliver of the plot and is hard to see.
    if xlim is None:
        xlim = xr

    # fit all requested models
    outs = {}
    for pf in prefactors:
        outs[pf] = fit_cavity_prefactor_compare(
            x, y, xr=xr, yc_um=yc_um, fit_yc=fit_yc,
            prefactor=pf, R_nm=R_nm,
            a0=a0, fit_a=(pf == 'powerlaw'),
            win=win, prom=prom, lam0_guess=lam0_guess, edges=edges
        )

    # print summary
    for pf, out in outs.items():
        p, d, met = out['params'], out['derived'], out['metrics']

        if pf == 'hankel':
            damp_txt = f"q_imag={p['q_imag_um^-1']:.3f} um^-1 ({d['q_imag_cm^-1']:.2e} cm^-1)"
        else:
            damp_txt = f"alpha_env={p['alpha_env_um^-1']:.3f} um^-1 ({d['alpha_env_cm^-1']:.2e} cm^-1)"

        print(f"{pf:8s}  q={d['q_cm^-1']:.2e} cm^-1  lambda={p['lambda_p_um']:.3f} um  "
              f"{damp_txt}  a={p['a']}  rmse={met['rmse']:.4g}  aic={met['aic']:.2f}")

    # sort by AIC
    order = sorted(outs.keys(), key=lambda k: outs[k]['metrics']['aic'])

    # plot
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=figsize, sharex=True,
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
    )

    ax.plot(x, y, 'k.', ms=4, alpha=0.7, label='data', zorder=5)

    # Fixed color per prefactor (not per AIC rank) so e.g. 'hankel' is always the
    # same color across every wavenumber's plot, regardless of which model wins AIC.
    prefactor_colors = {
        'hankel': '#1c7293', '1/sqrtx': '#e08214', 'none': '#7f7f7f',
        '1/x': '#2ca02c', 'powerlaw': '#9467bd',
    }

    for pf in order:
        out = outs[pf]
        xf, yf, m = out['x_fit_um'], out['y_fit'], out['mask']

        # Smooth dense curve for display; residuals below still use the real data grid.
        color = prefactor_colors.get(pf)
        line, = ax.plot(out['x_dense_um'], out['y_dense'], lw=1.6, color=color, label=pf)
        axr.plot(xf, y[m] - yf, lw=1.2, color=line.get_color())

    ax.set(
        title=f"{col}   (sorted by AIC)",
        ylabel='Amplitude (a.u.)',
        ylim=ylim
    )
    ax.set_ylabel('Amplitude (a.u.)', fontweight='bold')
    ax.legend(fontsize=12, ncol=2, frameon=False)
    ax.tick_params(direction='in', top=True, right=True)

    axr.axhline(0, color='gray', ls='--', lw=1)
    axr.set(
        xlim=xlim,
        xlabel=r'Distance ($\mu$m)',
        ylabel='resid.'
    )
    axr.set_xlabel(r'Distance ($\mu$m)', fontweight='bold')
    axr.tick_params(direction='in', top=True, right=True)

    if show_text:
        lines = []
        for pf in order:
            o = outs[pf]
            p = o['params']
            lam = p.get('lambda_p_um', 0)*1000
            q_re = 2*np.pi / (lam/1000) if lam > 0 else 0
            
            if pf == 'hankel':
                q_im = p.get('q_imag_um^-1', 1e-6)
            else:
                q_im = p.get('alpha_env_um^-1', 1e-6)
                
            damping = q_re / q_im if q_im > 1e-6 else np.inf
            rmse = o['metrics']['rmse']
            aic = o['metrics']['aic']
            
            lines.append(f"{pf}: λ={lam:.1f}nm, q={q_re:.2f}, γ⁻¹={damping:.1f}, RMSE={rmse:.3f}, AIC={aic:.1f}")
        
        txt = "\n".join(lines)
        ax.text(
            0.02, 0.02, txt, transform=ax.transAxes, fontsize=10,
            va='bottom', ha='left',

            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='0.8', alpha=0.9)
        )

    return outs, fig, (ax, axr)



def fit_two_hankel_shared_q(x_um, y1, y2, xr_fit=(0.2, 3.7), edges_um=(0.1, 3.7),
                            R_nm=25.0, share_q_imag=False, robust=True, win=15, prom=0.02, lam0_guess=None):
    """
    Joint fit of two channels with shared lambda_p (shared q), using a symmetric Hankel cavity model.

    Model per channel c:
        s_c(y) = B_c + A_c * Re{ exp(i phi_c) [H0^(1)(2 q_p,c xL) + H0^(1)(2 q_p,c xR)] }_norm

    with
        q_p,c = q + i*q_imag,c
        q = 2*pi/lambda_p

    Notes
    -----
    - The argument uses 2*q_p*x to represent a round-trip tip-edge-tip pathway.
    - Therefore the amplitude attenuation vs distance-from-edge x is asymptotically exp(-2*q_imag*x).
    - q is shared between O1A and O2A.
    - q_imag can be separate per channel, or shared if share_q_imag=True.
    """
    x_um = np.asarray(x_um, float)
    y1 = np.asarray(y1, float)
    y2 = np.asarray(y2, float)

    # fit mask
    m = (x_um >= xr_fit[0]) & (x_um <= xr_fit[1])
    x = x_um[m]
    s1 = y1[m]
    s2 = y2[m]

    # physical cavity edges
    left_edge, right_edge = edges_um
    R_um = R_nm * 1e-3
    eps = max(R_um, 1e-6)
    xL = np.clip(x - left_edge, eps, None)
    xR = np.clip(right_edge - x, eps, None)

    # initial lambda guess from standing-wave spacing: adjacent maxima ~ lambda/2
    def guess_lambda(xx, yy):
        wwin = win if (win % 2 == 1) else win + 1
        if len(yy) < max(wwin, 5):
            return None
        ys = savgol_filter(yy, wwin, 2)
        pks, _ = find_peaks(ys, prominence=prom)
        if len(pks) >= 2:
            dmax = np.median(np.diff(xx[pks]))
            return float(np.clip(2 * dmax, 0.05, 20.0))
        return None

    if lam0_guess is not None:
        lam0 = float(lam0_guess)
    else:
        lam_guesses = [guess_lambda(x, s1), guess_lambda(x, s2)]
        lam_guesses = [v for v in lam_guesses if v is not None]
        lam0 = np.median(lam_guesses) if lam_guesses else 3.0

    # rough per-channel init
    B10 = np.nanmedian(s1)
    A10 = 0.5 * (np.nanpercentile(s1, 95) - np.nanpercentile(s1, 5))
    B20 = np.nanmedian(s2)
    A20 = 0.5 * (np.nanpercentile(s2, 95) - np.nanpercentile(s2, 5))
    qimag10 = 0.5   # Im(q_p,1) in 1/um
    qimag20 = 0.5   # Im(q_p,2) in 1/um
    phi10, phi20 = 0.0, 0.0

    # parameterization
    # share_q_imag=False: [lam, B1, A1, qimag1, phi1, B2, A2, qimag2, phi2]
    # share_q_imag=True : [lam, qimag, B1, A1, phi1, B2, A2, phi2]
    if share_q_imag:
        p0 = np.array([lam0, 0.5, B10, A10, phi10, B20, A20, phi20], float)
        lb = np.array([1e-4, 0.0, -np.inf, -np.inf, -np.pi, -np.inf, -np.inf, -np.pi], float)
        ub = np.array([1e3, 50.0,  np.inf,  np.inf,  np.pi,  np.inf,  np.inf,  np.pi], float)
    else:
        p0 = np.array([lam0, B10, A10, qimag10, phi10, B20, A20, qimag20, phi20], float)
        lb = np.array([1e-4, -np.inf, -np.inf, 0.0, -np.pi, -np.inf, -np.inf, 0.0, -np.pi], float)
        ub = np.array([1e3,   np.inf,  np.inf, 50.0,  np.pi,  np.inf,  np.inf, 50.0,  np.pi], float)

    # scale residuals so one channel does not dominate
    s1_scale = max(np.std(s1), 1e-6)
    s2_scale = max(np.std(s2), 1e-6)

    def hankel_channel(lam, q_imag, B, A, phi):
        """q_p = q + i*q_imag, with 2*q_p*x for round-trip tip-launched path."""
        q = 2 * np.pi / lam
        q_p = q + 1j * q_imag
        h = hankel1(0, 2 * q_p * xL) + hankel1(0, 2 * q_p * xR)
        h = np.real(np.exp(1j * phi) * h)
        h = h / max(np.nanmax(np.abs(h)), 1e-12)   # normalize shape
        return B + A * h

    def unpack(theta):
        if share_q_imag:
            lam, q_imag, B1, A1, phi1, B2, A2, phi2 = theta
            qimag1 = qimag2 = q_imag
        else:
            lam, B1, A1, qimag1, phi1, B2, A2, qimag2, phi2 = theta
        return lam, (B1, A1, qimag1, phi1), (B2, A2, qimag2, phi2)

    def residual(theta):
        lam, c1, c2 = unpack(theta)
        f1 = hankel_channel(lam, c1[2], c1[0], c1[1], c1[3])
        f2 = hankel_channel(lam, c2[2], c2[0], c2[1], c2[3])
        r1 = (f1 - s1) / s1_scale
        r2 = (f2 - s2) / s2_scale
        return np.r_[r1, r2]

    res = least_squares(
        residual, p0, bounds=(lb, ub),
        loss='soft_l1' if robust else 'linear',
        f_scale=1.0, max_nfev=40000
    )

    # final outputs
    lam, c1, c2 = unpack(res.x)
    B1, A1, qimag1, phi1 = c1
    B2, A2, qimag2, phi2 = c2

    y1_fit = hankel_channel(lam, qimag1, B1, A1, phi1)
    y2_fit = hankel_channel(lam, qimag2, B2, A2, phi2)

    # raw residual metrics
    r1_raw = s1 - y1_fit
    r2_raw = s2 - y2_fit
    rss1 = np.sum(r1_raw**2)
    rss2 = np.sum(r2_raw**2)
    n1, n2 = len(s1), len(s2)
    rss_tot = rss1 + rss2
    n_tot = n1 + n2
    k = len(res.x)

    rmse1 = np.sqrt(rss1 / max(n1, 1))
    rmse2 = np.sqrt(rss2 / max(n2, 1))
    rmse_tot = np.sqrt(rss_tot / max(n_tot, 1))
    aic_tot = n_tot * np.log(max(rss_tot / n_tot, 1e-30)) + 2 * k

    q_um = 2 * np.pi / lam

    out = {
        'mask': m,
        'x_fit_um': x,
        'fit1': y1_fit,
        'fit2': y2_fit,
        'params_shared': {
            'lambda_p_um': lam,
            'q_rad_per_um': q_um,
            'q_cm^-1': q_um * 1e4,
            'share_q_imag': share_q_imag,
            'R_nm': R_nm,
            'xr_fit': xr_fit,
            'edges_um': edges_um
        },
        'params_ch1': {
            'B': B1, 'A': A1, 'q_imag_um^-1': qimag1, 'phi': phi1
        },
        'params_ch2': {
            'B': B2, 'A': A2, 'q_imag_um^-1': qimag2, 'phi': phi2
        },
        'derived': {
            'q_cm^-1': q_um * 1e4,
            'q_imag1_cm^-1': qimag1 * 1e4,
            'q_imag2_cm^-1': qimag2 * 1e4,
            # because amplitude ~ exp(-2*q_imag*x) in this round-trip convention:
            'amp_decay_len1_um': 1 / max(2 * qimag1, 1e-12),
            'amp_decay_len2_um': 1 / max(2 * qimag2, 1e-12),
        },
        'metrics': {
            'rmse1': rmse1, 'rmse2': rmse2, 'rmse_total': rmse_tot,
            'rss1': rss1, 'rss2': rss2, 'rss_total': rss_tot,
            'aic_total': aic_tot
        },
        'success': res.success,
        'fit_object': res
    }

    if share_q_imag:
        out['params_shared']['q_imag_um^-1'] = qimag1
        out['derived']['q_imag_cm^-1'] = qimag1 * 1e4
        out['derived']['amp_decay_len_um'] = 1 / max(2 * qimag1, 1e-12)

    return out


def fit_joint_hankel(amplp, w, xr_fit=(0.1, 3.7), edges_um=(0.1, 3.7),
                    R_nm=25.0, share_q_imag=True, win=7, prom=0.02,
                    xlim=(0.1, 3.7), figsize=(9.0, 5.2), ms=6, show_text=True, lam0_guess=None):
    """
    Joint Hankel fit for O1A and O2A at one frequency, with shared q.

    Parameters
    ----------
    amplp : DataFrame
        Must contain 'distance_um', '{w}cm-1_O1A', '{w}cm-1_O2A'
    w : str
        Example: '340.83'
    """
    col1 = f'{w}cm-1_O1A'
    col2 = f'{w}cm-1_O2A'

    x = amplp['distance_um'].values
    y1 = amplp[col1].values
    y2 = amplp[col2].values

    out_joint = fit_two_hankel_shared_q(
        x, y1, y2,
        xr_fit=xr_fit,
        edges_um=edges_um,
        R_nm=R_nm,
        share_q_imag=share_q_imag,
        win=win, prom=prom, lam0_guess=lam0_guess
    )

    ps = out_joint['params_shared']
    p1 = out_joint['params_ch1']
    p2 = out_joint['params_ch2']
    d = out_joint['derived']
    met = out_joint['metrics']

    print(f"shared q = {ps['q_cm^-1']:.2e} cm^-1   lambda = {ps['lambda_p_um']:.3f} um")
    if 'q_imag_um^-1' in ps:
        print(f"shared q_imag = {ps['q_imag_um^-1']:.3f} um^-1")
    print(f"O1A: q_imag={p1['q_imag_um^-1']:.3f} um^-1  rmse={met['rmse1']:.4g}")
    print(f"O2A: q_imag={p2['q_imag_um^-1']:.3f} um^-1  rmse={met['rmse2']:.4g}")
    print(f"total rmse={met['rmse_total']:.4g}  total AIC={met['aic_total']:.2f}")

    m = out_joint['mask']
    xf = out_joint['x_fit_um']
    f1 = out_joint['fit1']
    f2 = out_joint['fit2']

    fig, axs = plt.subplots(
        2, 2, figsize=figsize, sharex='col',
        gridspec_kw={'height_ratios':[3,1], 'hspace':0.05, 'wspace':0.18}
    )

    # top row
    axs[0,0].plot(x, y1, 'k.', ms=ms, alpha=0.7, label='O1A data', zorder=5)
    axs[0,0].plot(xf, f1, lw=1.6, label='joint hankel fit')
    axs[0,0].set(title=col1, ylabel='Amplitude (a.u.)', xlim=xlim)
    axs[0,0].legend(fontsize=8)

    axs[0,1].plot(x, y2, 'k.', ms=ms, alpha=0.7, label='O2A data', zorder=5)
    axs[0,1].plot(xf, f2, lw=1.6, label='joint hankel fit')
    axs[0,1].set(title=col2, xlim=xlim)
    axs[0,1].legend(fontsize=8)

    # bottom row
    axs[1,0].plot(xf, y1[m] - f1, lw=1.2)
    axs[1,0].axhline(0, color='gray', ls='--', lw=1)
    axs[1,0].set(xlabel=r'Distance ($\mu$m)', ylabel='resid.', xlim=xlim)

    axs[1,1].plot(xf, y2[m] - f2, lw=1.2)
    axs[1,1].axhline(0, color='gray', ls='--', lw=1)
    axs[1,1].set(xlabel=r'Distance ($\mu$m)', xlim=xlim)

    if show_text:
        if 'q_imag_um^-1' in ps:
            txt = (
                f"shared q = {ps['q_cm^-1']:.2e} cm$^{{-1}}$   "
                f"$\\lambda_p$ = {ps['lambda_p_um']:.3f} um\n"
                f"shared $q_{{imag}}$ = {ps['q_imag_um^-1']:.3f} um$^{{-1}}$\n"
                f"RMSE: O1A={met['rmse1']:.4g}, O2A={met['rmse2']:.4g}"
            )
        else:
            txt = (
                f"shared q = {ps['q_cm^-1']:.2e} cm$^{{-1}}$   "
                f"$\\lambda_p$ = {ps['lambda_p_um']:.3f} um\n"
                f"O1A $q_{{imag}}$={p1['q_imag_um^-1']:.3f} um$^{{-1}}$, "
                f"O2A $q_{{imag}}$={p2['q_imag_um^-1']:.3f} um$^{{-1}}$\n"
                f"RMSE: O1A={met['rmse1']:.4g}, O2A={met['rmse2']:.4g}"
            )

        axs[0,0].text(
            0.10, 0.02, txt, transform=axs[0,0].transAxes, fontsize=8,
            va='bottom', ha='left',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='0.8', alpha=0.9)
        )

    return out_joint, fig, axs

## Example ##
# out_joint, fig, axs = fit_joint_hankel(amplp, '344.17')
def plot_nf_fft(amplp, base_col, xr=(0.1, 3.7), 
                subtract_mean=True, fit_y=None, window='hann', pad_to=None,
                normalize_fft=True, plot_derivative=False,
                xlim_q=None, q_guess=None, figsize=(8, 4)):
    """
    Perform FFT analysis on a near-field line profile to extract polariton momentum (q).
    Modified to handle both O1 and O2 channels simultaneously based on the base label.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import get_window

    x_full = amplp['distance_um'].values
    
    # Identify channels (O1 and O2) based on base_col (e.g. '369.17' or '369.17cm-1')
    col_O1, col_O2 = None, None
    for suffix in ["_O1A", "_O1P"]:
        if f"{base_col}{suffix}" in amplp.columns:
            col_O1 = f"{base_col}{suffix}"
            col_O2 = f"{base_col}{suffix.replace('O1', 'O2')}"
            break
        elif f"{base_col}cm-1{suffix}" in amplp.columns:
            col_O1 = f"{base_col}cm-1{suffix}"
            col_O2 = f"{base_col}cm-1{suffix.replace('O1', 'O2')}"
            break
            
    if col_O1 is None or col_O2 not in amplp.columns:
        # Fallback if only one column or exactly as specified
        col_O1 = base_col
        col_O2 = base_col

    y1_full = amplp[col_O1].values
    y2_full = amplp[col_O2].values if col_O1 != col_O2 else y1_full

    # Step 1: Truncate data to xr
    mask = (x_full >= xr[0]) & (x_full <= xr[1])
    x = x_full[mask]
    dx = np.mean(np.diff(x))

    # Helper function to process single channel
    def _process_channel(y_full_c, fit_y_c):
        y_c = y_full_c[mask]
        if fit_y_c is not None:
            if len(fit_y_c) != len(y_c):
                raise ValueError(f"Length of fit_y ({len(fit_y_c)}) must match truncated data length ({len(y_c)}) within xr {xr}")
            y_proc_c = y_c - fit_y_c
            if subtract_mean:
                y_proc_c = y_proc_c - np.mean(y_proc_c)
        elif subtract_mean:
            y_proc_c = y_c - np.mean(y_c)
        else:
            y_proc_c = np.copy(y_c)
            
        def process_and_fft(data):
            if window:
                win_array = get_window(window, len(data))
                data_win = data * win_array
            else:
                data_win = data
                
            if pad_to is not None and pad_to > len(data_win):
                data_pad = np.zeros(pad_to)
                data_pad[:len(data_win)] = data_win
                n_fft_pts = pad_to
            else:
                data_pad = data_win
                n_fft_pts = len(data_win)
                
            fft_vals = np.fft.rfft(data_pad)
            amp_vals = np.abs(fft_vals)
            if normalize_fft and np.max(amp_vals) > 0:
                amp_vals = amp_vals / np.max(amp_vals)
            freqs = np.fft.rfftfreq(n_fft_pts, d=dx)
            qs = 2 * np.pi * freqs
            return data_win, qs, amp_vals

        y_win_c, q_um_c, fft_amp_c = process_and_fft(y_proc_c)
        y_deriv_c = np.gradient(y_proc_c, dx)
        y_deriv_win_c, _, fft_deriv_amp_c = process_and_fft(y_deriv_c)
        
        return y_c, y_proc_c, y_win_c, q_um_c, fft_amp_c, y_deriv_c, y_deriv_win_c, fft_deriv_amp_c

    # Handle fit_y which could be a tuple of two arrays or a single array
    fit_y1, fit_y2 = None, None
    if fit_y is not None:
        if isinstance(fit_y, tuple) and len(fit_y) == 2:
            fit_y1, fit_y2 = fit_y
        else:
            fit_y1 = fit_y
            fit_y2 = fit_y

    y1, y_proc1, y_win1, q_um, fft_amp1, y_deriv1, y_deriv_win1, fft_deriv_amp1 = _process_channel(y1_full, fit_y1)
    y2, y_proc2, y_win2, _, fft_amp2, y_deriv2, y_deriv_win2, fft_deriv_amp2 = _process_channel(y2_full, fit_y2)

    # Helper for peak finding
    def find_fft_peaks(amp_array, q_array, q_guess_c=None):
        if q_guess_c is not None:
            q_mask = (q_array >= q_guess_c * 0.8) & (q_array <= q_guess_c * 1.2)
        elif xlim_q:
            q_mask = (q_array >= xlim_q[0]) & (q_array >= 0.2) & (q_array <= xlim_q[1]) 
        else:
            q_mask = (q_array >= 0.5) & (q_array <= 30)
            
        if not np.any(q_mask): return []
        
        q_search = q_array[q_mask]
        amp_search = amp_array[q_mask]
        
        from scipy.signal import find_peaks
        pks, props = find_peaks(amp_search, prominence=np.max(amp_search)*0.05)
        
        if len(pks) == 0:
            return [q_search[np.argmax(amp_search)]]
            
        pks_sorted = pks[np.argsort(amp_search[pks])[::-1]]
        return q_search[pks_sorted].tolist()

    q_guess1, q_guess2 = q_guess, q_guess
    if isinstance(q_guess, (list, tuple, np.ndarray)) and len(q_guess) >= 2:
        q_guess1, q_guess2 = q_guess[0], q_guess[1]

    q_peaks1_orig = find_fft_peaks(fft_amp1, q_um, q_guess1)
    q_peaks1_deriv = find_fft_peaks(fft_deriv_amp1, q_um, q_guess1)
    q_peaks2_orig = find_fft_peaks(fft_amp2, q_um, q_guess2)
    q_peaks2_deriv = find_fft_peaks(fft_deriv_amp2, q_um, q_guess2)

    # Plotting
    n_cols = 3 if plot_derivative else 2
    n_rows = 2 if col_O1 != col_O2 else 1
    
    # Adjust figsize if using dual channels
    fig_w, fig_h = figsize
    if n_rows == 2 and fig_h <= 5:
        fig_h = fig_h * 2
    if n_cols == 3 and fig_w <= 8:
        fig_w = fig_w * 1.5
        
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), squeeze=False)

    def _plot_row(row_idx, col_name, y_raw, y_proc, y_win, fft_amp, y_deriv, y_deriv_win, fft_deriv_amp, pks_orig, pks_deriv, fit_y_c):
        label_orig = 'Original (mean sub)' if subtract_mean else 'Original Data'
        if fit_y_c is not None:
            axs[row_idx, 0].plot(x, y_raw, 'k.-', alpha=0.3, label='Raw Data')
            axs[row_idx, 0].plot(x, fit_y_c, 'orange', lw=1.5, ls='--', alpha=0.8, label='Fitted Curve')
            axs[row_idx, 0].plot(x, y_proc[:len(x)], 'm-', lw=1.5, label='Residuals (FFT Input)')
        else:
            axs[row_idx, 0].plot(x, y_proc, 'k.-', alpha=0.5, label=label_orig)
            if window:
                axs[row_idx, 0].plot(x, y_win[:len(x)], 'r-', lw=1.5, label=f'Windowed ({window})')
                
        title_suffix = '\n(Residual Analysis)' if fit_y_c is not None else ''
        axs[row_idx, 0].set(xlabel='Distance ($\mu$m)', ylabel='Amplitude (a.u.)', title=f'Spatial Signal {col_name}{title_suffix}')
        axs[row_idx, 0].legend(fontsize=8)
        
        # FFT Panel (position depends on if we plot derivative in middle)
        fft_idx = 2 if plot_derivative else 1
        axs[row_idx, fft_idx].plot(q_um, fft_amp, 'bo-', lw=1.5, ms=4)
        axs[row_idx, fft_idx].set(xlabel='q ($\mu$m$^{-1}$)', ylabel='FFT Amp' + (' (norm)' if normalize_fft else ''), title='Momentum Spectrum')
        if xlim_q: axs[row_idx, fft_idx].set_xlim(xlim_q)
        else: axs[row_idx, fft_idx].set_xlim(0, 30)

        if len(pks_orig) > 0:
            q_max = pks_orig[0]
            axs[row_idx, fft_idx].axvline(q_max, color='red', linestyle='--', alpha=0.6)
            idx_max = np.argmin(np.abs(q_um - q_max))
            axs[row_idx, fft_idx].text(q_max*1.05, fft_amp[idx_max]*0.9, f'q={q_max:.2f}', color='red', fontsize=9)
            
        if plot_derivative:
            axs[row_idx, 1].plot(x, y_deriv, 'g.-', alpha=0.5, label='Derivative (dy/dx)')
            if window:
                axs[row_idx, 1].plot(x, y_deriv_win[:len(x)], 'r-', lw=1.5, label=f'Windowed')
            axs[row_idx, 1].set(xlabel='Distance ($\mu$m)', ylabel='d(Amp)/dx', title='Spatial Derivative')
            axs[row_idx, 1].legend(fontsize=8)
            
    _plot_row(0, col_O1, y1, y_proc1, y_win1, fft_amp1, y_deriv1, y_deriv_win1, fft_deriv_amp1, q_peaks1_orig, q_peaks1_deriv, fit_y1)
    if n_rows == 2:
        _plot_row(1, col_O2, y2, y_proc2, y_win2, fft_amp2, y_deriv2, y_deriv_win2, fft_deriv_amp2, q_peaks2_orig, q_peaks2_deriv, fit_y2)

    plt.tight_layout()
    
    wavenumber_str = base_col.split('cm-1')[0] if 'cm-1' in base_col else base_col
    
    out_dict = {
        'x_um': x,
        'q_um^-1': q_um,
        f'{wavenumber_str}_fft': {
            'qfft_O1': q_peaks1_orig[0] if q_peaks1_orig else None,
            'qfft_O2': q_peaks2_orig[0] if q_peaks2_orig else None,
            'qfft_deriv_O1': q_peaks1_deriv[0] if q_peaks1_deriv else None,
            'qfft_deriv_O2': q_peaks2_deriv[0] if q_peaks2_deriv else None,
            # 保留一个 peakvalue 以防其他旧代码使用，但去掉后面那些冗长的列表
            'peakvalue': q_peaks1_orig
        }
    }
        
    return fig, axs, out_dict


def append_nf_results(results_df, wn, outs_joint=None, out_dict_fft=None, outs_cavity=None, out_joint=None):
    """
    Helper function to extract fitted parameters and FFT peaks, appending them to a pandas DataFrame.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        The DataFrame to append to. If None, a new one is created.
    wn : float or str
        The wavenumber for this row.
    outs_joint : dict or None
        A dictionary of output dictionaries from joint_fit, e.g., {'hankel': out_hankel, '1/sqrtx': out_sqrtx}.
    out_dict_fft : dict or None
        The output dictionary from plot_nf_fft.
    outs_cavity : dict or None
        The output dictionary from compare_cavity_models.
    out_joint : dict or None
        (Legacy) single output dictionary. Stored as 'hankel' generic.
        
    Returns
    -------
    pd.DataFrame
        The updated DataFrame containing the new row.
    """
    import pandas as pd
    
    # Initialize DataFrame if None
    if results_df is None or results_df.empty:
        results_df = pd.DataFrame()
        
    if 'wavenumber (cm^-1)' not in results_df.columns:
        results_df['wavenumber (cm^-1)'] = pd.Series(dtype=float)
        
    try:
        wn_val = float(wn)
    except ValueError:
        wn_val = wn
        
    new_data = {'wavenumber (cm^-1)': wn_val}
    
    # 1. Extract from Fit Dictionaries
    # Handle single out_joint legacy
    if out_joint is not None and outs_joint is None:
        outs_joint = {'hankel': out_joint}
        
    if outs_joint is not None:
        for model_name, out_j in outs_joint.items():
            try:
                lam_fit = q_fit = None
                if 'params_shared' in out_j:
                    lam_fit = out_j['params_shared']['lambda_p_um']
                    q_fit = out_j['params_shared']['q_rad_per_um']
                elif 'lambda_p_um' in out_j:
                    lam_fit = out_j['lambda_p_um']
                    q_fit = out_j['q_rad_per_um']
                    
                if lam_fit is not None:
                    new_data[f'lambda_joint_{model_name} (um)'] = lam_fit
                    new_data[f'q_joint_{model_name} (um^-1)'] = q_fit
            except Exception as e:
                print(f"Warning: Failed to extract fit parameters for {model_name}: {e}")

    # 2. Extract from FFT Dictionary
    if out_dict_fft is not None:
        try:
            try:
                wn_float_str = str(float(wn))
            except ValueError:
                wn_float_str = str(wn)
                
            matched_key = None
            for key in out_dict_fft.keys():
                if '_fft' in key and '_deriv_fft' not in key and key.replace('_fft', '').startswith(wn_float_str.rstrip('0').rstrip('.')):
                    matched_key = key
                    break
                    
            if not matched_key:
                fft_key = f"{str(wn).split('cm-1')[0]}_fft"
                if fft_key in out_dict_fft:
                    matched_key = fft_key

            if matched_key:
                fft_dict = out_dict_fft[matched_key]
                
                # New logic for dual-channel format
                if 'qfft_O1' in fft_dict and fft_dict['qfft_O1'] is not None:
                    new_data['q_fft_O1 (um^-1)'] = fft_dict['qfft_O1']
                if 'qfft_O2' in fft_dict and fft_dict['qfft_O2'] is not None:
                    new_data['q_fft_O2 (um^-1)'] = fft_dict['qfft_O2']
                if 'qfft_deriv_O1' in fft_dict and fft_dict['qfft_deriv_O1'] is not None:
                    new_data['q_deriv_fft_O1 (um^-1)'] = fft_dict['qfft_deriv_O1']
                if 'qfft_deriv_O2' in fft_dict and fft_dict['qfft_deriv_O2'] is not None:
                    new_data['q_deriv_fft_O2 (um^-1)'] = fft_dict['qfft_deriv_O2']
                
                # Legacy logic
                if 'qfft_O1' not in fft_dict:
                    peaks = fft_dict.get('peakvalue', [])
                    if len(peaks) > 0:
                        new_data['q_fft (um^-1)'] = peaks[0]
                        
                    deriv_key = matched_key.replace('_fft', '_deriv_fft')
                    if deriv_key in out_dict_fft:
                        deriv_peaks = out_dict_fft[deriv_key].get('peakvalue', [])
                        if len(deriv_peaks) > 0:
                            new_data['q_deriv_fft (um^-1)'] = deriv_peaks[0]
                        
        except Exception as e:
            print(f"Warning: Failed to extract FFT parameters: {e}")

    # 3. Extract from compare_cavity_models (outs_cavity)
    if outs_cavity is not None:
        for method, model_out in outs_cavity.items():
            try:
                if 'derived' in model_out and 'q_rad_per_um' in model_out['derived']:
                    q_val = model_out['derived']['q_rad_per_um']
                    col_name = f'q_{method}_fit (um^-1)'
                    new_data[col_name] = q_val
            except Exception as e:
                print(f"Warning: Failed to extract cavity model {method}: {e}")

    # Update or Append avoiding duplicates
    mask = results_df['wavenumber (cm^-1)'] == wn_val
    if mask.any():
        idx = results_df.index[mask].tolist()[0]
        for col_name, val in new_data.items():
            if val is not None:
                results_df.loc[idx, col_name] = val
    else:
        new_row_df = pd.DataFrame([new_data])
        if not new_row_df.isna().all().all():
            results_df = pd.concat([results_df, new_row_df], ignore_index=True)
            
    cols = ['wavenumber (cm^-1)'] + [c for c in results_df.columns if c != 'wavenumber (cm^-1)']
    return results_df[cols]


def joint_fit(amplp, w, model='hankel',
              edges='double',                          # MANDATORY: 'single' or 'double'
              xr_fit=(0.1, 3.7), edges_um=(0.1, 3.7),
              R_nm=25.0, share_q_imag=True, win=7, prom=0.02,
              robust=True, lam0_guess=None,
              xlim=(0.1, 3.7), figsize=(9.0, 5.2), ms=6, show_text=True):
    import numpy as np
    from scipy.optimize import least_squares
    from scipy.signal import savgol_filter, find_peaks
    from scipy.special import hankel1
    import matplotlib.pyplot as plt
 
    if edges not in ('single', 'double'):
        raise ValueError("edges must be 'single' or 'double'")
 
    col1 = f'{w}cm-1_O1A'
    col2 = f'{w}cm-1_O2A'
 
    x  = amplp['distance_um'].values
    y1 = amplp[col1].values
    y2 = amplp[col2].values
 
    R_um = R_nm * 1e-3
    eps  = max(R_um, 1e-6)
 
    def _edge_distances(xf):
        xL = np.clip(xf - edges_um[0], eps, None)
        xR = np.clip(edges_um[1] - xf, eps, None) if edges == 'double' else None
        return xL, xR
 
    def _geom_prefactor(xL, xR):
        if edges == 'single':
            G = 1.0 / np.sqrt(xL + R_um)
        else:
            G = 1.0/np.sqrt(xL + R_um) + 1.0/np.sqrt(xR + R_um)
        return G / max(np.nanmax(np.abs(G)), 1e-12)
 
    def _envelope(alpha, xf, xL):
        if edges == 'single':
            env = np.exp(-alpha * xL)
        else:
            yc  = 0.5 * (edges_um[0] + edges_um[1])
            env = np.cosh(alpha * (xf - yc))
        return env / max(np.nanmax(np.abs(env)), 1e-12)
 
    def _lam_guess(xf, s):
        ys_g   = savgol_filter(s, 5, 2)
        pks, _ = find_peaks(ys_g, prominence=0.005)
        dmax0  = np.median(np.diff(xf[pks])) if len(pks) >= 2 else 0.4
        lam0   = float(np.clip(2*dmax0, 0.05, 20.0))
        return lam0, 0.7*lam0, 1.2*lam0
 
    # ================================================================== #
    #  HANKEL                                                             #
    # ================================================================== #
    if model == 'hankel':
 
        if edges == 'double':
            out = fit_two_hankel_shared_q(
                x, y1, y2,
                xr_fit=xr_fit, edges_um=edges_um,
                R_nm=R_nm, share_q_imag=share_q_imag,
                win=win, prom=prom, lam0_guess=lam0_guess
            )
            ps  = out['params_shared']
            p1  = out['params_ch1']
            p2  = out['params_ch2']
            met = out['metrics']
            m   = out['mask']
            xf  = out['x_fit_um']
            f1  = out['fit1']
            f2  = out['fit2']
 
            print(f"shared q = {ps['q_cm^-1']:.2e} cm^-1   lambda = {ps['lambda_p_um']:.3f} um")
            if 'q_imag_um^-1' in ps:
                print(f"shared q_imag = {ps['q_imag_um^-1']:.3f} um^-1")
            print(f"O1A: q_imag={p1['q_imag_um^-1']:.3f} um^-1  rmse={met['rmse1']:.4g}")
            print(f"O2A: q_imag={p2['q_imag_um^-1']:.3f} um^-1  rmse={met['rmse2']:.4g}")
            print(f"total rmse={met['rmse_total']:.4g}  total AIC={met['aic_total']:.2f}")
 
            if show_text:
                if 'q_imag_um^-1' in ps:
                    txt = (
                        f"shared q={ps['q_cm^-1']:.2e} cm$^{{-1}}$, "
                        f"$\\lambda_p$={ps['lambda_p_um']:.3f} um\n"
                        f"shared $q_{{imag}}$={ps['q_imag_um^-1']:.3f} um$^{{-1}}$\n"
                        f"RMSE: O1A={met['rmse1']:.4g}, O2A={met['rmse2']:.4g}"
                    )
                else:
                    txt = (
                        f"shared q={ps['q_cm^-1']:.2e} cm$^{{-1}}$, "
                        f"$\\lambda_p$={ps['lambda_p_um']:.3f} um\n"
                        f"O1A $q_{{imag}}$={p1['q_imag_um^-1']:.3f}, "
                        f"O2A $q_{{imag}}$={p2['q_imag_um^-1']:.3f} um$^{{-1}}$\n"
                        f"RMSE: O1A={met['rmse1']:.4g}, O2A={met['rmse2']:.4g}"
                    )
 
        else:
            m  = (x >= xr_fit[0]) & (x <= xr_fit[1])
            xf = x[m]
            s1 = y1[m];  s2 = y2[m]
            xL, _ = _edge_distances(xf)
 
            lam0, lam_lo, lam_hi = _lam_guess(xf, s1)
            if lam0_guess is not None:
                lam0 = float(lam0_guess)
                lam_lo = min(lam_lo, lam0 * 0.7)
                lam_hi = max(lam_hi, lam0 * 1.3)
 
            def _hankel_single(xL, B, A, q_imag, lam, phi):
                q   = 2*np.pi / lam
                q_p = q + 1j*q_imag
                h   = hankel1(0, 2*q_p*xL)
                hre = np.real(np.exp(1j*phi) * h)
                hre = hre / max(np.nanmax(np.abs(hre)), 1e-12)
                return B + A * hre
 
            def _ch_p0(s):
                return [np.nanmedian(s),
                        0.5*(np.nanpercentile(s,95) - np.nanpercentile(s,5)),
                        0.5, 0.0]
 
            p1_0 = _ch_p0(s1);  p2_0 = _ch_p0(s2)
            p0 = np.array([lam0, 0.5,
                           p1_0[0], p1_0[1], p1_0[3],
                           p2_0[0], p2_0[1], p2_0[3]])
            lb = np.array([lam_lo, 0.0,
                           -np.inf, -np.inf, -np.pi,
                           -np.inf, -np.inf, -np.pi])
            ub = np.array([lam_hi, 50.0,
                            np.inf,  np.inf,  np.pi,
                            np.inf,  np.inf,  np.pi])
 
            def _resid_h1(theta):
                lam, qi        = theta[0], theta[1]
                B1, A1, phi1   = theta[2:5]
                B2, A2, phi2   = theta[5:8]
                r1 = _hankel_single(xL, B1, A1, qi, lam, phi1) - s1
                r2 = _hankel_single(xL, B2, A2, qi, lam, phi2) - s2
                return np.concatenate([r1, r2])
 
            res = least_squares(
                _resid_h1, p0, bounds=(lb, ub),
                loss='soft_l1' if robust else 'linear',
                f_scale=max(np.std(s1)*0.3, 1e-6),
                max_nfev=30_000
            )
 
            lam_fit, qi_fit  = res.x[0], res.x[1]
            B1, A1, phi1     = res.x[2:5]
            B2, A2, phi2     = res.x[5:8]
 
            f1 = _hankel_single(xL, B1, A1, qi_fit, lam_fit, phi1)
            f2 = _hankel_single(xL, B2, A2, qi_fit, lam_fit, phi2)
 
            q_um = 2*np.pi / lam_fit
            n, k = 2*len(xf), len(res.x)
 
            def _met(fit, s):
                rss  = np.sum((fit-s)**2)
                rmse = np.sqrt(rss / max(len(s)-k, 1))
                aic  = len(s)*np.log(max(rss/len(s), 1e-30)) + 2*k
                return dict(rss=rss, rmse=rmse, aic=aic)
 
            met1     = _met(f1, s1);  met2 = _met(f2, s2)
            rss_tot  = met1['rss'] + met2['rss']
            rmse_tot = np.sqrt(rss_tot / max(n-k, 1))
            aic_tot  = n*np.log(max(rss_tot/n, 1e-30)) + 2*k
 
            print(f"shared q = {q_um*1e4:.2e} cm^-1   lambda = {lam_fit:.3f} um")
            print(f"shared q_imag = {qi_fit:.3f} um^-1")
            print(f"O1A: rmse={met1['rmse']:.4g}")
            print(f"O2A: rmse={met2['rmse']:.4g}")
            print(f"total rmse={rmse_tot:.4g}  total AIC={aic_tot:.2f}")
 
            out = dict(
                model         = 'hankel_single',
                mask          = m, x_fit_um=xf, fit1=f1, fit2=f2,
                params_shared = {'lambda_p_um': lam_fit, 'q_rad_per_um': q_um, 'q_cm^-1': q_um*1e4,
                                 'q_imag_um^-1': qi_fit, 'R_nm': R_nm},
                params_ch1    = {'B': B1, 'A': A1, 'phi': phi1},
                params_ch2    = {'B': B2, 'A': A2, 'phi': phi2},
                metrics       = {'rmse1': met1['rmse'], 'rmse2': met2['rmse'],
                                 'rmse_total': rmse_tot, 'aic_total': aic_tot},
                success       = res.success, fit_object=res
            )
            ps  = out['params_shared']
            met = out['metrics']
 
            if show_text:
                txt = (
                    f"shared q={ps['q_cm^-1']:.2e} cm$^{{-1}}$, "
                    f"$\\lambda_p$={ps['lambda_p_um']:.3f} um\n"
                    f"shared $q_{{imag}}$={ps['q_imag_um^-1']:.3f} um$^{{-1}}$\n"
                    f"RMSE: O1A={met['rmse1']:.4g}, O2A={met['rmse2']:.4g}"
                )
 
        fit_label = 'hankel fit'
 
    # ================================================================== #
    #  1/sqrtx                                                            #
    # ================================================================== #
    elif model == '1/sqrtx':
        m  = (x >= xr_fit[0]) & (x <= xr_fit[1])
        xf = x[m]
        s1 = y1[m];  s2 = y2[m]
 
        xL, xR = _edge_distances(xf)
        G       = _geom_prefactor(xL, xR)
 
        lam0, lam_lo, lam_hi = _lam_guess(xf, s1)
        if lam0_guess is not None:
            lam0 = float(lam0_guess)
            lam_lo = min(lam_lo, lam0 * 0.7)
            lam_hi = max(lam_hi, lam0 * 1.3)
 
        def _sqrtx_model(xf, G, B, A, alpha, lam, phi):
            q   = 2*np.pi / lam
            env = _envelope(alpha, xf, xL)
            return B + A * np.sin(2*q*xf + phi) * env * G
 
        def _ch_p0(s):
            return [np.nanmedian(s),
                    0.5*(np.nanpercentile(s,95) - np.nanpercentile(s,5)),
                    0.5, 0.0]
 
        p1_0 = _ch_p0(s1);  p2_0 = _ch_p0(s2)
        p0   = np.array([lam0] + p1_0 + p2_0)
        lb   = np.array([lam_lo, -np.inf, -np.inf, 0.0, -np.pi,
                                  -np.inf, -np.inf, 0.0, -np.pi])
        ub   = np.array([lam_hi,  np.inf,  np.inf, 50., np.pi,
                                   np.inf,  np.inf, 50., np.pi])
 
        def _resid(theta):
            lam              = theta[0]
            B1, A1, a1, phi1 = theta[1:5]
            B2, A2, a2, phi2 = theta[5:9]
            r1 = _sqrtx_model(xf, G, B1, A1, a1, lam, phi1) - s1
            r2 = _sqrtx_model(xf, G, B2, A2, a2, lam, phi2) - s2
            return np.concatenate([r1, r2])
 
        res = least_squares(
            _resid, p0, bounds=(lb, ub),
            loss='soft_l1' if robust else 'linear',
            f_scale=max(np.std(s1)*0.3, 1e-6),
            max_nfev=30_000
        )
 
        lam_fit          = res.x[0]
        B1,A1,a1,phi1    = res.x[1:5]
        B2,A2,a2,phi2    = res.x[5:9]
 
        f1 = _sqrtx_model(xf, G, B1, A1, a1, lam_fit, phi1)
        f2 = _sqrtx_model(xf, G, B2, A2, a2, lam_fit, phi2)
 
        q_um = 2*np.pi / lam_fit
        n, k = 2*len(xf), len(res.x)
 
        def _met(fit, s):
            rss  = np.sum((fit-s)**2)
            rmse = np.sqrt(rss / max(len(s)-k, 1))
            aic  = len(s)*np.log(max(rss/len(s), 1e-30)) + 2*k
            return dict(rss=rss, rmse=rmse, aic=aic)
 
        met1     = _met(f1, s1);  met2 = _met(f2, s2)
        rss_tot  = met1['rss'] + met2['rss']
        rmse_tot = np.sqrt(rss_tot / max(n-k, 1))
        aic_tot  = n*np.log(max(rss_tot/n, 1e-30)) + 2*k
 
        print(f"shared lambda = {lam_fit:.3f} um   q = {q_um*1e4:.2e} cm^-1")
        print(f"O1A: alpha_env={a1:.3f} um^-1  rmse={met1['rmse']:.4g}")
        print(f"O2A: alpha_env={a2:.3f} um^-1  rmse={met2['rmse']:.4g}")
        print(f"total rmse={rmse_tot:.4g}  total AIC={aic_tot:.2f}")
 
        out = dict(
            model         = '1/sqrtx',
            mask          = m, x_fit_um=xf, fit1=f1, fit2=f2, G=G,
            params_shared = {'lambda_p_um': lam_fit, 'q_rad_per_um': q_um, 'q_cm^-1': q_um*1e4, 'R_nm': R_nm},
            params_ch1    = {'B': B1, 'A': A1, 'alpha_env_um^-1': a1, 'phi': phi1},
            params_ch2    = {'B': B2, 'A': A2, 'alpha_env_um^-1': a2, 'phi': phi2},
            metrics       = {'rmse1': met1['rmse'], 'rmse2': met2['rmse'],
                             'rmse_total': rmse_tot, 'aic_total': aic_tot},
            success       = res.success, fit_object=res
        )
        met = out['metrics']
        ps  = out['params_shared']
 
        fit_label = r'$1/\sqrt{x}$ fit'
 
        if show_text:
            txt = (
                f"shared $\\lambda_p$={lam_fit:.3f} um, "
                f"q={q_um*1e4:.2e} cm$^{{-1}}$\n"
                f"O1A $\\alpha_{{env}}$={a1:.3f} um$^{{-1}}$, "
                f"O2A $\\alpha_{{env}}$={a2:.3f} um$^{{-1}}$\n"
                f"RMSE: O1A={met1['rmse']:.4g}, O2A={met2['rmse']:.4g}"
            )
 
    else:
        raise ValueError(f"Unknown model '{model}'. Choose 'hankel' or '1/sqrtx'.")
 
    # ================================================================== #
    #  Shared plotting                                                    #
    # ================================================================== #
    m  = out['mask']
    xf = out['x_fit_um']
    f1 = out['fit1']
    f2 = out['fit2']
 
    fig, axs = plt.subplots(
        2, 2, figsize=figsize, sharex='col',
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05, 'wspace': 0.18}
    )
 
    axs[0,0].plot(x,  y1, 'k.', ms=ms, alpha=0.7, label='O1A data', zorder=5)
    axs[0,0].plot(xf, f1, lw=1.6, label=fit_label)
    axs[0,0].set(title=col1, ylabel='Amplitude (a.u.)', xlim=xlim)
    axs[0,0].legend(fontsize=8)
 
    axs[0,1].plot(x,  y2, 'k.', ms=ms, alpha=0.7, label='O2A data', zorder=5)
    axs[0,1].plot(xf, f2, lw=1.6, label=fit_label)
    axs[0,1].set(title=col2, xlim=xlim)
    axs[0,1].legend(fontsize=8)
 
    axs[1,0].plot(xf, y1[m] - f1, lw=1.2)
    axs[1,0].axhline(0, color='gray', ls='--', lw=1)
    axs[1,0].set(xlabel=r'Distance ($\mu$m)', ylabel='resid.', xlim=xlim)
 
    axs[1,1].plot(xf, y2[m] - f2, lw=1.2)
    axs[1,1].axhline(0, color='gray', ls='--', lw=1)
    axs[1,1].set(xlabel=r'Distance ($\mu$m)', xlim=xlim)
 
    if show_text:
        axs[0,0].text(
            0.10, 0.02, txt, transform=axs[0,0].transAxes, fontsize=8,
            va='bottom', ha='left',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='0.8', alpha=0.9)
        )
 
    return out, fig, axs



# ==========================================
# 🌌 Near-Field Exclusive FFT Dispersion Analysis Engine (Pure Version)
# ==========================================
def plot_nf_fft_2d(imgdata, label, ch, 
                   fft_mode='amp',        # 'amp', 'phase', or 'complex'
                   xr=None,               # Spatial range for FFT extraction [min_um, max_um]
                   mirror_at=False,       # False, float, or 'auto' to enable spatial folding and averaging 
                   probe_freqs=None,      # Frequencies (list) used for 'auto' correlation (e.g., [334, 340])
                   preview_fold=None,     # Show 2D folded map? None means auto-show only when fft_mode == 'amp'
                   ax_preview=None,       # Provide (ax_amp, ax_phase) to plot preview on existing axes
                   clim_fft=None,         # FFT color limit [vmin, vmax], None for auto
                   xlim_k=[-2.5, 2.5],    # X-axis (Wavevector q) display range [min, max], None for auto
                   ylim_freq=None,        # Y-axis (Frequency) display range [min, max], None for auto
                   padding_factor=1,      # Zero-padding factor for FFT (enhances k-axis resolution)
                   window_type='hann',    # Window function type: 'hann' or 'boxcar'
                   pad_mode='zero',       # 'zero' for zero-padding, 'mean' for padding with edge means
                   remove_dc=True,        # Whether to explicitly remove the mean of each spatial slice
                   apply_1st_deriv=False, # Apply first derivative filter along q-axis
                   apply_2nd_deriv=False, # Apply second derivative filter
                   deriv_sigma=2.0,       # Gaussian smoothing radius for second derivative
                   savefig=False,
                   registry=None):        # Added registry parameter to access save paths
    """
    Receives imgdata and plots a highly customizable FFT dispersion map.
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.fft import fft, fftfreq, fftshift
    from scipy.signal import windows
    from scipy.ndimage import gaussian_filter1d

    # ==========================================
    # 🎨 Customization Section
    # ==========================================
    # -- 1. Canvas and Axes Control --
    fig_size = (4, 6)         
    axes_ticks_dir = 'in'     
    
    # 🌟 Core: Axis ranges and units
    freq_unit = 'cm-1'      # Swap units: 'cm-1' or 'THz' (converts both X and Y automatically)
    
    # -- 2. Title and Color Contrast --
    # 🌟 Core: Modify title here
    fold_str = f" (Folded\u2728)" if (mirror_at is not False and mirror_at is not None) else ""
    custom_title = f"{ch} {fft_mode.capitalize()} FFT{fold_str}" if ch in ["O1", "O2"] else f"{ch} FFT{fold_str}" 
    custom_position = None  # Title position [X, Y], None for auto center
    
    # cmap_fft = 'magma'    # Recommended: 'magma', 'inferno', 'viridis', etc.
    cmap_fft = 'magma'
    
    # -- 3. Colorbar Global Control --
    cbar_orientation = 'horizontal'       # 'vertical' or 'horizontal'
    cbar_rect = [0.165, 0.89, 0.7, 0.018] # Absolute positioning [X, Y, Width, Height]. None for auto
    cbar_outline_color = 'black'          # Border color, None to hide
    cbar_outline_width = 1.0      
    
    cbar_ticks_dir = 'in'           
    cbar_ticks_pos = 'top'          # 'top', 'bottom', 'left', 'right'
    cbar_label_pos = 'top'          
    
    cbar_font_size = 10           
    cbar_title_fontsize = 12      
    cbar_title_pad = -10            # Negative pulls title towards colorbar
    
    cbar_shrink = 0.8             
    cbar_aspect = 15            
    cbar_pad = 0.05               
    # ==========================================

    # 1. Unpack data smoothly accommodating different dictionary keys natively from plotnanoFTIR
    if 'freq_cm' in imgdata:
        freq_cm = imgdata['freq_cm']
    elif 'wavenumber' in imgdata:
        freq_cm = imgdata['wavenumber'][:, 0] if imgdata['wavenumber'].ndim == 2 else imgdata['wavenumber']
    else:
        raise ValueError("Frequency array not found in imgdata!")

    if 'dist_um' in imgdata:
        dist = imgdata['dist_um']
    elif 'distance_um' in imgdata:
        dist = imgdata['distance_um'][0, :] if imgdata['distance_um'].ndim == 2 else imgdata['distance_um']
    else:
        raise ValueError("Distance array not found in imgdata!")
        
    amp_source = imgdata.get('amp', imgdata.get(f'{ch}A'))
    phase_source = imgdata.get('phase', imgdata.get(f'{ch}P'))
    
    if amp_source is None or phase_source is None:
        raise ValueError(f"Amplitude/Phase data not found for channel {ch}.")
        
    amp_vals = amp_source.values if hasattr(amp_source, 'values') else amp_source
    phase_vals = phase_source.values if hasattr(phase_source, 'values') else phase_source
    
    if fft_mode == 'amp':
        target_vals = amp_vals
    elif fft_mode == 'phase':
        target_vals = phase_vals
    elif fft_mode == 'complex':
        target_vals = amp_vals * np.exp(1j * phase_vals)
    else:
        raise ValueError("fft_mode must be 'amp', 'phase', or 'complex'")

    # 2. Select spatial ROI
    mask = np.ones(len(dist), dtype=bool)
    if xr is not None:
        mask = (dist >= xr[0]) & (dist <= xr[1])
    dist_roi = dist[mask]
    target_roi = target_vals[:, mask]
    
    if len(dist_roi) < 2:
        print("Error: Too few data points in ROI for FFT!")
        return imgdata, None

    # --- 🌟 Fold Spatial ROI if Requested ---
    if mirror_at is not False and mirror_at is not None:
        N_r = len(dist_roi)
        nf = 0
        xc_idx = 0
        if mirror_at == 'auto':
            if probe_freqs is None:
                probe_idx = np.arange(len(freq_cm))
            else:
                probe_idx = [np.argmin(np.abs(freq_cm - f)) for f in probe_freqs]
                
            best_score = -np.inf
            best_xc_idx = N_r // 2
            best_nf = 1
            best_l, best_r = 0, N_r - 1
            
            # We let the algorithm trim up to 25% of the garbage points on either edge to maximize cavity symmetry 
            max_l = max(1, N_r // 4)
            min_r = min(N_r - 1, N_r - N_r // 4)
            
            for l_idx in range(max_l + 1):
                for r_idx in range(min_r, N_r):
                    if r_idx - l_idx < 10: 
                        continue
                    
                    c_idx = (l_idx + r_idx) // 2
                    cur_nf = min(c_idx - l_idx, r_idx - c_idx)
                    if cur_nf < 5: 
                        continue
                    
                    total_r = 0.0
                    n_v = 0
                    for ti in probe_idx:
                        row = target_roi[ti, :]
                        if np.iscomplexobj(row):
                            row_mag = np.abs(row)
                            lv = row_mag[c_idx - cur_nf : c_idx][::-1]
                            rv = row_mag[c_idx + 1 : c_idx + 1 + cur_nf]
                        else:
                            lv = row[c_idx - cur_nf : c_idx][::-1]
                            rv = row[c_idx + 1 : c_idx + 1 + cur_nf]
                        
                        lv_ = lv - np.mean(lv)
                        rv_ = rv - np.mean(rv)
                        d = np.std(lv_) * np.std(rv_)
                        if d > 1e-12:
                            total_r += np.dot(lv_, rv_) / (cur_nf * d)
                            n_v += 1
                    
                    score = total_r / n_v if n_v > 0 else np.nan
                    if not np.isnan(score) and score > best_score:
                        best_score = score
                        best_xc_idx = c_idx
                        best_nf = cur_nf
                        best_l = l_idx
                        best_r = r_idx
                        
            xc_idx = best_xc_idx
            nf = best_nf
            print(f"[Auto Fold 2D Smart-Crop] Adjusted range: {dist_roi[best_l]:.2f} to {dist_roi[best_r]:.2f} µm, Center: {dist_roi[xc_idx]:.3f} µm (Score: {best_score:.3f}, Pts: {nf})")
        else:
            try:
                xc = float(mirror_at)
                xc_idx = np.argmin(np.abs(dist_roi - xc))
                left_len = xc_idx
                right_len = len(dist_roi) - xc_idx - 1
                nf = min(left_len, right_len)
                if nf < 1:
                    print(f"⚠️ Warning: mirror_at={xc} is too close to the edge! Skipping fold.")
                else:
                    print(f"[Manual Fold] at {dist_roi[xc_idx]:.3f} µm, Pts: {nf}")
            except ValueError:
                print(f"⚠️ Warning: mirror_at={mirror_at} is not valid. Skipping fold.")
                nf = 0

        if nf >= 1:
            amp_target = amp_vals[:, mask]
            phase_target = phase_vals[:, mask]
            
            folded_amp = np.zeros((len(freq_cm), nf), dtype=float)
            folded_phase = np.zeros((len(freq_cm), nf), dtype=float)
            
            for i in range(len(freq_cm)):
                # Amp
                lv_a = amp_target[i, xc_idx - nf : xc_idx][::-1]
                rv_a = amp_target[i, xc_idx + 1 : xc_idx + 1 + nf]
                folded_amp[i] = 0.5 * (lv_a + rv_a)
                # Phase
                lv_p = phase_target[i, xc_idx - nf : xc_idx][::-1]
                rv_p = phase_target[i, xc_idx + 1 : xc_idx + 1 + nf]
                folded_phase[i] = 0.5 * (lv_p + rv_p)
                
            dist_roi = dist_roi[xc_idx + 1 : xc_idx + 1 + nf] - dist_roi[xc_idx]
            
            if fft_mode == 'amp':
                target_roi = folded_amp
            elif fft_mode == 'phase':
                target_roi = folded_phase
            elif fft_mode == 'complex':
                target_roi = folded_amp * np.exp(1j * folded_phase)
            
            # --- 🌟 Preview the Folded 2D Spatial Map ---
            if preview_fold or (preview_fold is None and fft_mode == 'amp'):
                if ax_preview is not None:
                    ax1, ax2 = ax_preview
                    fig_sp = ax1.figure
                else:
                    fig_sp, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 4), sharey=True)
                    
                ff, dd = np.meshgrid(freq_cm, dist_roi)
                
                try:
                    from snippet import Sky
                    cmap_a = Sky
                except ImportError:
                    cmap_a = 'viridis'
                    
                im1 = ax1.pcolormesh(ff, dd, folded_amp.T, cmap=cmap_a, shading='auto')
                fig_sp.colorbar(im1, ax=ax1, label='Amplitude')
                ax1.set_title(f"Folded {ch}A")
                ax1.set_xlabel('Freq. (cm$^{-1}$)')
                ax1.set_ylabel('Dist. (\u03bcm)')
                
                im2 = ax2.pcolormesh(ff, dd, folded_phase.T, cmap='RdBu', shading='auto')
                fig_sp.colorbar(im2, ax=ax2, label='Phase')
                ax2.set_title(f"Folded {ch}P")
                ax2.set_xlabel('Freq. (cm$^{-1}$)')
                if ax_preview is None:
                    ax2.set_ylabel('Dist. (\u03bcm)')
                
                if ax_preview is None:
                    plt.tight_layout()
                    plt.show()  # Display inline without disrupting return values

    # 3. FFT Calculation
    # We must match the Old Code's spatial orientation: edge on left, center on right.
    # folded_amp currently has center on left (index 0). We reverse it.
    if mirror_at is not None:
        target_roi = target_roi[:, ::-1]

    N = target_roi.shape[1]
    N_pad = int(N * padding_factor) if N > 0 else 1024 
    dx = np.mean(np.diff(dist_roi))
    k_vals = 2 * np.pi * fftshift(fftfreq(N_pad, d=dx))
    
    spatial_bg_profile = np.mean(target_roi, axis=0)
    
    fft_map = []
    fft_map_complex = []
    if window_type == 'hann':
        win = windows.hann(N)
    else:
        win = windows.boxcar(N)
    
    for i in range(len(freq_cm)):
        E = target_roi[i, :]
        #E = E - spatial_bg_profile 
        if remove_dc:
            E = E - np.mean(E) 
            
        E_win = E * win
        
        if N_pad > N:
            if pad_mode == 'mean':
                n_avg = min(5, N)
                left_mean = np.mean(E_win[:n_avg])
                right_mean = np.mean(E_win[-n_avg:])
                pad_left = (N_pad - N) // 2
                pad_right = N_pad - N - pad_left
                E_padded = np.pad(E_win, (pad_left, pad_right), mode='constant', constant_values=(left_mean, right_mean))
            elif pad_mode == 'edge':
                pad_left = (N_pad - N) // 2
                pad_right = N_pad - N - pad_left
                E_padded = np.pad(E_win, (pad_left, pad_right), mode='edge')
            else:
                # pad_mode == 'zero'
                # Old code specifically appends physical 0.0, not the tail mean.
                E_padded = np.pad(E_win, (0, N_pad - N), mode='constant', constant_values=0.0)
        else:
            E_padded = E_win
            
        ft = fftshift(fft(E_padded, n=N_pad))
        fft_map.append(np.abs(ft))
        fft_map_complex.append(ft)
        
    fft_map = np.array(fft_map)
    fft_map_complex = np.array(fft_map_complex)
    
    # 💥 Keep k=0 so the dot falls exactly on 0, but manually flatten its height to 0.0 to prevent scaling issues.
    mask_pos = (k_vals >= 0)
    k_vals = k_vals[mask_pos]
    
    fft_map = fft_map[:, mask_pos]
    fft_map_complex = fft_map_complex[:, mask_pos]
    
    # Set the DC literal value to 0 so it plots at q=0 but with 0 intensity
    if len(k_vals) > 0 and k_vals[0] == 0:
        fft_map[:, 0] = 0.0

    # 4. Derivative processing
    fft_plot_data = fft_map.copy()
    default_cbar_label = 'Intensity (a.u.)'
    _npad_str = f"N={N}→{N_pad}(×{padding_factor})" if N_pad > N else f"N={N}"
    default_title = f'{ch} Dispersion ({fft_mode.capitalize()})  [{_npad_str}]'


    if apply_1st_deriv:
        d1_fft = gaussian_filter1d(fft_map, sigma=deriv_sigma, axis=1, order=1)
        fft_plot_data = d1_fft
        # Force symmetric colorbar and exclude k=0 jump artifacts
        c_idx = d1_fft.shape[1] // 2
        exc = max(1, int(d1_fft.shape[1] * 0.02)) # exclude center 2%
        valid_d1 = np.concatenate((d1_fft[:, :c_idx-exc], d1_fft[:, c_idx+exc:]), axis=1)
        max_val = np.max(np.abs(valid_d1)) if valid_d1.size > 0 else np.max(np.abs(d1_fft))
        clim_fft = [-max_val, max_val]
        cmap_fft = 'RdBu_r'
        default_cbar_label = 'dI/dq'
        default_title = f'{ch} Dispersion (1st Deriv q)'
    elif apply_2nd_deriv:
        d2_fft = gaussian_filter1d(fft_map, sigma=deriv_sigma, axis=0, order=2)
        fft_plot_data = -d2_fft
        fft_plot_data[fft_plot_data < 0] = 0 
        default_cbar_label = '-d²I/dE²'
        default_title = f'{ch} Dispersion (2nd Deriv)'

    # 5. 🌟 Physical units bidirectional conversion core 🌟
    to_eV_from_cm = 1 / 8065.544
    to_THz_from_cm = 0.0299792458
    to_THz_from_um_inv = 299.792458 # speed of light c ≈ 299.79 um/ps
    
    if freq_unit == 'THz':
        # Y-axis conversion
        y_vals = freq_cm * to_THz_from_cm
        ylabel = 'Frequency (THz)'
        sec_label = 'Energy (meV)'
        func_forward = lambda x: x * 4.135667
        func_inverse = lambda x: x / 4.135667
        
        # X-axis conversion
        x_vals = k_vals * to_THz_from_um_inv
        xlabel = r'Wavevector $qc$ (THz)'
        # Smart adaptation for display range
        xlim_plot = [x * to_THz_from_um_inv for x in xlim_k] if xlim_k else None
    else:
        # Keep Y-axis
        y_vals = freq_cm
        ylabel = r'Frequency (cm$^{-1}$)'
        sec_label = 'Energy (meV)'
        func_forward = lambda x: x * to_eV_from_cm * 1000
        func_inverse = lambda x: x / (to_eV_from_cm * 1000)
        
        # Keep X-axis
        x_vals = k_vals
        xlabel = r'Wavevector $q$ ($\rm{\mu m}^{-1}$)'
        xlim_plot = xlim_k

    # 6. Plotting
    fig_fft, ax_fft = plt.subplots(figsize=fig_size)
    
    vmin = clim_fft[0] if clim_fft else np.min(fft_plot_data)
    vmax = clim_fft[1] if clim_fft else np.max(fft_plot_data)
    
    kk, yy = np.meshgrid(x_vals, y_vals)
    # Changed rasterised parameter for broader compatibility
    im = ax_fft.pcolormesh(kk, yy, fft_plot_data, cmap=cmap_fft, shading='auto', 
                           vmin=vmin, vmax=vmax, rasterized=True)
    
    ax_fft.set_xlabel(xlabel)
    ax_fft.set_ylabel(ylabel)
    ax_fft.set_title(custom_title if custom_title else default_title, y=1.08, loc='center' if custom_position is None else None)
    if custom_position is not None:
        ax_fft.title.set_position(custom_position)
    
    if xlim_plot: ax_fft.set_xlim(xlim_plot)
    if ylim_freq: ax_fft.set_ylim(ylim_freq)
    ax_fft.tick_params(axis='both', direction=axes_ticks_dir)

    secax_fft = ax_fft.secondary_yaxis('right', functions=(func_forward, func_inverse))
    secax_fft.set_ylabel(sec_label, color='gray')
    secax_fft.tick_params(axis='y', direction=axes_ticks_dir, colors='gray')

    # 7. Colorbar
    if cbar_rect is not None:
        cax = fig_fft.add_axes(cbar_rect)
        cb = fig_fft.colorbar(im, cax=cax, orientation=cbar_orientation)
    else:
        cb = fig_fft.colorbar(im, ax=ax_fft, orientation=cbar_orientation,
                              shrink=cbar_shrink, aspect=cbar_aspect, pad=cbar_pad)
    
    cb.set_ticks([vmin, vmax])
    fmt = lambda x: f"{x:.2e}" if (abs(x) < 0.01 or abs(x) > 1000) and x != 0 else f"{x:.1f}"
    cb.set_ticklabels([fmt(vmin), fmt(vmax)])
    
    if cbar_outline_color:
        cb.outline.set_visible(True)
        cb.outline.set_edgecolor(cbar_outline_color)
        cb.outline.set_linewidth(cbar_outline_width)
    else:
        cb.outline.set_visible(False)
        
    cb.ax.tick_params(labelsize=cbar_font_size, direction=cbar_ticks_dir)
    if cbar_orientation == 'horizontal':
        cb.ax.xaxis.set_ticks_position(cbar_ticks_pos)
        cb.ax.xaxis.set_label_position(cbar_label_pos)
    else:
        cb.ax.yaxis.set_ticks_position(cbar_ticks_pos)
        cb.ax.yaxis.set_label_position(cbar_label_pos)
        
    cb.set_label(default_cbar_label, fontsize=cbar_title_fontsize, labelpad=cbar_title_pad)

    # 8. Save Data
    imgdata[f'fft_k_vals_{fft_mode}'] = k_vals  # Save raw k (for downstream alignment)
    imgdata[f'fft_map_{fft_mode}'] = fft_plot_data
    
    if savefig and registry is not None and label in registry:
        try:
            prefix = registry[label]["NearField"]["channels"][ch]["save_path_prefix"]
            save_dir = os.path.dirname(prefix)
            if save_dir:
                os.makedirs(save_dir, exist_ok=True)
            suffix = "2ndDeriv" if apply_2nd_deriv else "Raw"
            out_path = f"{prefix}_FFT_{fft_mode}_{suffix}.png"
            fig_fft.savefig(out_path, dpi=600, bbox_inches='tight')
            print(f"✅ Saved to: {out_path}")
        except KeyError:
            print("⚠️ Warning: registry does not contain expected keys for save path.")

    return imgdata, fig_fft


def plot_stacked_fft_from_2d(imgdata_O1, imgdata_O2,
                              target_freqs,
                              # ── same FFT params as plot_nf_fft_2d ──
                              xr=None,
                              mirror_at=False,
                              probe_freqs=None,
                              padding_factor=2,
                              window_type='boxcar',
                              remove_dc=True,
                              # ── display params ──
                              q_range=(0, 30),
                              q_guess=None,
                              search_window=0.2,
                              stacked_figsize=(11, 7.5),
                              stacked_offset=0.3,
                              save_dir=None,
                              label='Sample',
                              save_peaks=True,
                              manual_peaks_dict=None,
                              power_scale=1.0):
    """
    Compute Amp / Phase / Complex FFT for a set of target wavenumbers
    from the raw 2D imgdata, and display them as a 3-panel stacked waterfall.

    Uses the **same** fold, padding, and windowing logic as plot_nf_fft_2d
    (just pass the same parameters), so the slices are perfectly consistent
    with the 2D dispersion map — no drift at high frequencies.

    Parameters
    ----------
    imgdata_O1/O2 : dict
        Raw imgdata from plotnanoFTIR (or from df['O1']/df['O2']).
        Must contain amplitude/phase 2D arrays and distance/frequency axes.
    target_freqs : list of float
        Wavenumber values (cm-1) to slice.
    xr, mirror_at, probe_freqs, padding_factor, window_type, remove_dc :
        Identical meaning as in plot_nf_fft_2d.
    q_range : (qmin, qmax) µm-1
    q_guess : list, length == len(target_freqs).
        Each element: scalar, list, or None (skip peak annotation).
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from scipy.fft import fft, fftfreq, fftshift
    from scipy.signal import windows as sig_windows
    import os

    def _compute_fft_slices(imgdata, ch):
        """Run fold+window+FFT on the full 2D data, return slices at target_freqs."""
        # ── 1. Unpack ──
        if 'freq_cm' in imgdata:
            freq_cm = np.array(imgdata['freq_cm'])
        elif 'wavenumber' in imgdata:
            wn = imgdata['wavenumber']
            freq_cm = wn[:, 0] if (hasattr(wn, 'ndim') and wn.ndim == 2) else np.array(wn)
        else:
            raise KeyError(f"Frequency axis not found in imgdata for {ch}")

        if 'distance_um' in imgdata:
            dist = imgdata['distance_um']
            dist = dist[0, :] if (hasattr(dist, 'ndim') and dist.ndim == 2) else np.array(dist)
        else:
            raise KeyError(f"Distance axis not found in imgdata for {ch}")

        amp_src   = imgdata.get('amp',   imgdata.get(f'{ch}A'))
        phase_src = imgdata.get('phase', imgdata.get(f'{ch}P'))
        if amp_src is None or phase_src is None:
            raise ValueError(f"Amp/Phase not found for channel {ch}.")

        amp_v   = amp_src.values   if hasattr(amp_src,   'values') else np.array(amp_src)
        phase_v = phase_src.values if hasattr(phase_src, 'values') else np.array(phase_src)

        # ── 2. Spatial ROI ──
        mask_d = np.ones(len(dist), dtype=bool)
        if xr is not None:
            mask_d = (dist >= xr[0]) & (dist <= xr[1])
        dist_roi   = dist[mask_d]
        amp_roi    = amp_v[:,   mask_d]
        phase_roi  = phase_v[:, mask_d]

        # ── 3. Folding (mirrors plot_nf_fft_2d exactly) ──
        nf_pts = 0
        if mirror_at is not False and mirror_at is not None:
            N_r = len(dist_roi)
            if mirror_at == 'auto':
                pidx = ([np.argmin(np.abs(freq_cm - f)) for f in probe_freqs]
                        if probe_freqs else np.arange(len(freq_cm)))
                best_score, best_xc, best_nf = -np.inf, N_r // 2, 1
                max_l = max(1, N_r // 4)
                min_r = min(N_r - 1, N_r - N_r // 4)
                for l_i in range(max_l + 1):
                    for r_i in range(min_r, N_r):
                        if r_i - l_i < 10:
                            continue
                        c_i  = (l_i + r_i) // 2
                        cnf  = min(c_i - l_i, r_i - c_i)
                        if cnf < 5:
                            continue
                        tot, n_v = 0.0, 0
                        for ti in pidx:
                            row = amp_roi[ti, :]
                            lv = row[c_i - cnf : c_i][::-1] ; rv = row[c_i + 1 : c_i + 1 + cnf]
                            lv_ = lv - np.mean(lv) ; rv_ = rv - np.mean(rv)
                            d = np.std(lv_) * np.std(rv_)
                            if d > 1e-12:
                                tot += np.dot(lv_, rv_) / (cnf * d) ; n_v += 1
                        score = tot / n_v if n_v > 0 else np.nan
                        if not np.isnan(score) and score > best_score:
                            best_score, best_xc, best_nf = score, c_i, cnf
                xc_idx, nf_pts = best_xc, best_nf
                print(f"[{ch} Auto-fold] center={dist_roi[xc_idx]:.3f} µm, nf={nf_pts}")
            else:
                xc_idx = np.argmin(np.abs(dist_roi - float(mirror_at)))
                nf_pts = min(xc_idx, len(dist_roi) - xc_idx - 1)
                print(f"[{ch} Manual fold] center={dist_roi[xc_idx]:.3f} µm, nf={nf_pts}")

            if nf_pts >= 1:
                folded_amp   = np.zeros((len(freq_cm), nf_pts))
                folded_phase = np.zeros((len(freq_cm), nf_pts))
                for i in range(len(freq_cm)):
                    lv_a = amp_roi[i,   xc_idx - nf_pts : xc_idx][::-1]
                    rv_a = amp_roi[i,   xc_idx + 1      : xc_idx + 1 + nf_pts]
                    folded_amp[i]   = 0.5 * (lv_a + rv_a)
                    lv_p = phase_roi[i, xc_idx - nf_pts : xc_idx][::-1]
                    rv_p = phase_roi[i, xc_idx + 1      : xc_idx + 1 + nf_pts]
                    folded_phase[i] = 0.5 * (lv_p + rv_p)
                dist_roi  = dist_roi[xc_idx + 1 : xc_idx + 1 + nf_pts]
                amp_roi   = folded_amp
                phase_roi = folded_phase

        # Reverse to match plot_nf_fft_2d convention (edge-left, center-right)
        if nf_pts >= 1:
            amp_roi   = amp_roi[:,   ::-1]
            phase_roi = phase_roi[:, ::-1]

        # ── 4. FFT ──
        N   = amp_roi.shape[1]
        N_p = max(int(N * padding_factor), N)
        dx  = float(np.mean(np.diff(dist_roi))) if len(dist_roi) > 1 else 1.0
        k_all = 2 * np.pi * fftshift(fftfreq(N_p, d=dx))

        win = sig_windows.hann(N) if window_type == 'hann' else sig_windows.boxcar(N)

        maps = {}   # mode -> 2D array shape (n_freq, N_p//2+1 or similar)
        for mode in ('amp', 'phase', 'complex'):
            if mode == 'amp':
                raw = amp_roi
            elif mode == 'phase':
                raw = phase_roi
            else:
                raw = amp_roi * np.exp(1j * phase_roi)

            fft_rows = []
            for i in range(len(freq_cm)):
                E = raw[i, :]
                if remove_dc:
                    E = E - np.mean(E)
                E_win = E * win
                E_pad = np.pad(E_win, (0, N_p - N), mode='constant', constant_values=0.0)
                ft = fftshift(fft(E_pad, n=N_p))
                fft_rows.append(np.abs(ft))
            maps[mode] = np.array(fft_rows)

        # Keep positive k only
        mask_pos = k_all >= 0
        k_pos    = k_all[mask_pos]
        for mode in maps:
            maps[mode] = maps[mode][:, mask_pos]
        if len(k_pos) > 0 and k_pos[0] == 0:
            for mode in maps:
                maps[mode][:, 0] = 0.0

        # ── 5. Extract target freq slices ──
        extracted = {}   # freq -> {actual_freq, q, amp_spec, phase_spec, complex_spec}
        q_mask = (k_pos > 0) & (k_pos <= q_range[1] * 10.0)  # q_range is in 10⁵ cm⁻¹, k_pos is in µm⁻¹
        q_trim = k_pos[q_mask]
        for f in target_freqs:
            row_idx   = int(np.argmin(np.abs(freq_cm - f)))
            actual_f  = float(freq_cm[row_idx])
            extracted[f] = {
                'actual_freq':   actual_f,
                'q':             q_trim,
                'amp_spec':      maps['amp'][row_idx,     q_mask],
                'phase_spec':    maps['phase'][row_idx,   q_mask],
                'complex_spec':  maps['complex'][row_idx, q_mask],
            }
        return extracted, N, N_p

    # ── Compute for each channel ──
    _r_O1 = _compute_fft_slices(imgdata_O1, 'O1') if imgdata_O1 is not None else ({}, 0, 0)
    _r_O2 = _compute_fft_slices(imgdata_O2, 'O2') if imgdata_O2 is not None else ({}, 0, 0)
    ext_O1, _N_O1, _Np_O1 = _r_O1
    ext_O2, _N_O2, _Np_O2 = _r_O2

    def _find_peak(q, spec, qg):
        lo, hi = qg * (1 - search_window), qg * (1 + search_window)
        m = (q >= lo) & (q <= hi)
        if not np.any(m):
            return None, None
        return float(q[m][np.argmax(spec[m])]), float(np.max(spec[m]))

    def _draw_stack(extracted, ch_label, N_orig=0, N_pad_val=0):
        sorted_f = sorted([f for f in target_freqs if f in extracted])
        if not sorted_f:
            return []

        _npad_tag = f"N={N_orig}→{N_pad_val}(×{padding_factor})" if N_pad_val > N_orig else f"N={N_orig}"
        fig, (ax_a, ax_p, ax_c) = plt.subplots(1, 3, figsize=stacked_figsize)
        fig.suptitle(f"{ch_label} stacked 1-D FFT  |  from 2D map  |  mirror={mirror_at}  [{_npad_tag}]",
                     fontsize=10, y=0.98)
        ax_a.set_title(f"{ch_label} Amplitude", pad=10)
        ax_p.set_title(f"{ch_label} Phase",     pad=10)
        ax_c.set_title(f"{ch_label} Complex FFT",    pad=10)

        cmap = cm.get_cmap('plasma')
        wn_min, wn_max = sorted_f[0], sorted_f[-1]
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=wn_min, vmax=wn_max))
        sm.set_array([])
        cbar_ax = fig.add_axes([0.92, 0.35, 0.015, 0.3])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Frequency (cm$^{-1}$)')

        peak_rows = []

        for idx, wn in enumerate(sorted_f):
            dat     = extracted[wn]
            q       = dat['q']
            q_plot  = q / 10.0            # µm⁻¹ → 10⁵ cm⁻¹ (1 µm⁻¹ = 10⁴ cm⁻¹ = 0.1×10⁵ cm⁻¹)
            sa      = dat['amp_spec']
            sp      = dat['phase_spec']
            sc      = dat['complex_spec']
            actual_f = dat['actual_freq']
            c        = cmap((wn - wn_min) / max(wn_max - wn_min, 1e-5))
            y_base   = idx * stacked_offset

            # Normalise independently (guard zero) and apply power scaling to enhance weak peaks
            def _nrm(s):
                mx = np.max(s)
                return (s / mx)**power_scale if mx > 0 else s

            ax_a.plot(q_plot, _nrm(sa) + y_base, '.-', ms=4, color=c, lw=1.5)
            ax_p.plot(q_plot, _nrm(sp) + y_base, '.-', ms=4, color=c, lw=1.5)
            ax_c.plot(q_plot, _nrm(sc) + y_base, '.-', ms=4, color=c, lw=1.5)

            # ── peak annotation ──
            found = []
            
            if manual_peaks_dict and wn in manual_peaks_dict:
                # Use predefined explicit peaks for this wavenumber
                key = f'Cmpx{ch_label}'
                if key in manual_peaks_dict[wn]:
                    explicit_qs = manual_peaks_dict[wn][key]
                    
                    # distinct colors for 1st, 2nd, 3rd, 4th... peak
                    peak_colors = ['#55a868', '#c44e52', '#8172b3', '#ccb974', '#4c72b0', '#64b5cd']
                    
                    for pi, p_q in enumerate(explicit_qs):
                        p_color = peak_colors[pi % len(peak_colors)]
                        idx_q = np.argmin(np.abs(q_plot - p_q))
                        ax_c.plot(p_q, _nrm(sc)[idx_q] + y_base, 'o', mec=p_color, mfc='none', mew=1.5, ms=7)
                        
                        # Vertical dashed lines removed as requested
                        
                        txt_va = 'bottom' if idx % 2 == 0 else 'top'
                        txt_dy = +0.04 if idx % 2 == 0 else -0.04
                        ax_c.text(p_q, _nrm(sc)[idx_q] + y_base + txt_dy, f"{p_q:.2f}",
                                  color=p_color, fontsize=8, ha='center', va=txt_va)
            else:
                # Use q_guess to search
                cur_guess = None
                if q_guess is not None and len(q_guess) == len(target_freqs):
                    elem = q_guess[target_freqs.index(wn)]
                    if elem is not None:
                        cur_guess = list(elem) if (hasattr(elem, '__iter__') and not isinstance(elem, str)) else [elem]

                if cur_guess:
                    for qg in cur_guess:
                        # search on amp spectrum (qg in µm⁻¹ → convert to 10⁵ cm⁻¹)
                        p_q, p_v = _find_peak(q_plot, _nrm(sa), qg / 10.0)
                        if p_q is not None:
                            found.append(p_q * 10.0)   # store back in µm⁻¹ for CSV
                            ax_a.plot(p_q, p_v + y_base, 'o', mec='#4c72b0', mfc='none', mew=1.5, ms=7)
                            # Alternate text above/below to reduce overlap
                            txt_va = 'bottom' if idx % 2 == 0 else 'top'
                            txt_dy = +0.04 if idx % 2 == 0 else -0.04
                            ax_a.text(p_q, p_v + y_base + txt_dy, f"{p_q:.2f}",
                                      color='#4c72b0', fontsize=8, ha='center', va=txt_va)
                        # same q on phase panel
                        if p_q is not None:
                            idx_q = np.argmin(np.abs(q_plot - p_q))
                            ax_p.plot(p_q, _nrm(sp)[idx_q] + y_base, 'o', mec='#dd8452', mfc='none', mew=1.5, ms=7)
                            ax_c.plot(p_q, _nrm(sc)[idx_q] + y_base, 'o', mec='#55a868', mfc='none', mew=1.5, ms=7)

            # freq label (placed at the right end of the data line, but inside the plot)
            for ax_i in (ax_a, ax_p, ax_c):
                # Put the label at the last data point, but aligned right so it extends leftwards over the line
                label_x = q_plot[-1]
                # Slightly above the line
                ax_i.text(label_x, y_base + 0.05, f"{actual_f:.0f}",
                          color=c, fontsize=8, fontweight='bold', va='bottom', ha='right')

            row = {'wn_target': wn, 'actual_freq': actual_f}
            for pi, pq in enumerate(found):
                row[f'peak_{pi+1}_q'] = pq
            peak_rows.append(row)

        q_xlim = (q_range[0], q_range[1])     # directly use 10⁵ cm⁻¹
        k_xlim = (q_xlim[0] * 10 / (2 * np.pi),
                  q_xlim[1] * 10 / (2 * np.pi))              # k = q/(2π) in µm⁻¹

        for i_ax, ax in enumerate((ax_a, ax_p, ax_c)):
            ax.set_xlim(q_xlim)
            ax.set_yticks([])
            ax.set_xlabel(r'$q = 2\pi/r$  ($10^5$ cm$^{-1}$)')
            ax.tick_params(direction='in', which='both', top=True, right=True)

            # twiny top axis: k = 1/r in µm⁻¹
            at = ax.twiny()
            at.set_xlim(k_xlim)
            at.set_xlabel(r'$k = 1/r$  ($\mu$m$^{-1}$)', fontsize=9)
            at.tick_params(direction='in', which='both')

        # y-axis label on leftmost panel only
        ax_a.set_ylabel('FFT amplitude (a.u.)', fontsize=9)

        plt.subplots_adjust(wspace=0.08)
        # Avoid double render in %matplotlib widget (ipympl) mode
        try:
            import matplotlib
            _bk = matplotlib.get_backend().lower()
            if 'widget' not in _bk and 'ipympl' not in _bk:
                try:
                    from IPython.display import display as _ipy_display
                    _ipy_display(fig)
                except ImportError:
                    plt.show()
        except Exception:
            plt.show()
        return peak_rows

    rows_O1 = _draw_stack(ext_O1, 'O1', _N_O1, _Np_O1)
    rows_O2 = _draw_stack(ext_O2, 'O2', _N_O2, _Np_O2)

    if save_dir and save_peaks:
        import pandas as pd
        os.makedirs(save_dir, exist_ok=True)
        for ch_label, rows in [('O1', rows_O1), ('O2', rows_O2)]:
            if rows:
                df_out = pd.DataFrame(rows)
                fname  = f"{label}_2DFFT_peaks_{ch_label}_pad{padding_factor}_m{str(mirror_at).replace('.','p')}.csv"
                fpath  = os.path.join(save_dir, fname)
                df_out.to_csv(fpath, index=False)
                print(f"✅ Saved → {fpath}")

    return {'O1': ext_O1, 'O2': ext_O2}

# ==========================================
# 🌟 Stacked Channel FFT System (Added for parallel O1/O2 analysis)
# ==========================================

from scipy.signal import get_window
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator

def plot_channel_fft(x, amp, phase, label, wn=None,
                     xr=(0.24, 3.4), 
                     mirror_at=None,
                     q_range=(0, 5),              # 10⁵ cm⁻¹, applies to all 3 FFT panels
                     window='boxcar', 
                     pad_factor=2.0,            
                     pad_mode='end',              
                     q_guess=[4, 22], 
                     search_window=0.3,
                     peak_method='find_peak',
                     tick_dir='in',    
                     y_bins=5,
                     plot=True):        
    
    # 提取全量数据
    x_full = x.values if isinstance(x, pd.Series) else np.array(x)
    amp_full = amp.values if isinstance(amp, pd.Series) else np.array(amp)
    phase_full = phase.values if isinstance(phase, pd.Series) else np.array(phase)

    # 截断数据 (真正用于FFT分析的数据)
    mask = (x_full >= xr[0]) & (x_full <= xr[1])
    x_c = x_full[mask]
    amp_c = amp_full[mask]
    phase_c = phase_full[mask]
    
    # --- 🌟 1D Folding Logic ---
    if mirror_at is not False and mirror_at is not None:
        N_r = len(x_c)
        if mirror_at == 'auto':
            best_score = -np.inf
            best_xc_idx = N_r // 2
            best_nf = 1
            best_l, best_r = 0, N_r - 1
            max_l = max(1, N_r // 4)
            min_r = min(N_r - 1, N_r - N_r // 4)
            for l_idx in range(max_l + 1):
                for r_idx in range(min_r, N_r):
                    if r_idx - l_idx < 10: continue
                    c_idx = (l_idx + r_idx) // 2
                    cur_nf = min(c_idx - l_idx, r_idx - c_idx)
                    if cur_nf < 5: continue
                    
                    lv_a = amp_c[c_idx - cur_nf : c_idx][::-1]
                    rv_a = amp_c[c_idx + 1 : c_idx + 1 + cur_nf]
                    lv_a_ = lv_a - np.mean(lv_a)
                    rv_a_ = rv_a - np.mean(rv_a)
                    da = np.std(lv_a_) * np.std(rv_a_)
                    score_a = np.dot(lv_a_, rv_a_) / (cur_nf * da) if da > 1e-12 else np.nan
                    
                    lv_p = phase_c[c_idx - cur_nf : c_idx][::-1]
                    rv_p = phase_c[c_idx + 1 : c_idx + 1 + cur_nf]
                    lv_p_ = lv_p - np.mean(lv_p)
                    rv_p_ = rv_p - np.mean(rv_p)
                    dp = np.std(lv_p_) * np.std(rv_p_)
                    score_p = np.dot(lv_p_, rv_p_) / (cur_nf * dp) if dp > 1e-12 else np.nan
                    
                    if not np.isnan(score_a) and not np.isnan(score_p):
                        score = (score_a + score_p) / 2
                        if score > best_score:
                            best_score, best_xc_idx, best_nf, best_l, best_r = score, c_idx, cur_nf, l_idx, r_idx
                            
            xc_idx, nf = best_xc_idx, best_nf
            print(f"[{label} Auto Fold] Cut: {x_c[best_l]:.2f} to {x_c[best_r]:.2f}, Center: {x_c[xc_idx]:.3f} µm (Pts: {nf})")
        else:
            xc = float(mirror_at)
            xc_idx = np.argmin(np.abs(x_c - xc))
            nf = min(xc_idx, len(x_c) - xc_idx - 1)
            
        if nf >= 1:
            lv_a = amp_c[xc_idx - nf : xc_idx][::-1]
            rv_a = amp_c[xc_idx + 1 : xc_idx + 1 + nf]
            amp_c = 0.5 * (lv_a + rv_a)[::-1]  # reverse to make center on right
            
            lv_p = phase_c[xc_idx - nf : xc_idx][::-1]
            rv_p = phase_c[xc_idx + 1 : xc_idx + 1 + nf]
            phase_c = 0.5 * (lv_p + rv_p)[::-1]
            
            x_c_unfolded = x_c[xc_idx + 1 : xc_idx + 1 + nf]
            # x=0 对应 Edge（arr[0]），x 增大方向 → Mirror（arr[-1]），保留原始间距
            x_c_plot = x_c_unfolded - x_c_unfolded[0]
            x_c = x_c_plot
            
            is_folded = True
    else:
        is_folded = False
    # ----------------------------

    dx = np.mean(np.diff(x_c)) if len(x_c) > 1 else 1.0
    
    # Amp、Phase 及 Complex 构建与去均值 (去背景)
    amp_proc = amp_c - np.mean(amp_c)
    phase_proc = phase_c - np.mean(phase_c)
    complex_data = amp_c * np.exp(1j * phase_c)
    complex_proc = complex_data - np.mean(complex_data)
    
    # 获取Window
    n_points = len(amp_proc)
    win_array = get_window(window, n_points) if window and window != 'boxcar' else np.ones(n_points)
    
    amp_win = amp_proc * win_array
    phase_win = phase_proc * win_array
    complex_win = complex_proc * win_array
    
    # Zero Padding 
    pad_to = int(np.ceil(n_points * pad_factor)) if pad_factor is not None and pad_factor >= 1.0 else n_points
    n_zeros = pad_to - n_points
    
    # 平滑 Padding 基准值
    n_avg = min(5, n_points)
    if n_points > 0:
        amp_end, phase_end, cplx_end = np.mean(amp_win[-n_avg:]), np.mean(phase_win[-n_avg:]), np.mean(complex_win[-n_avg:])
        amp_start, phase_start, cplx_start = np.mean(amp_win[:n_avg]), np.mean(phase_win[:n_avg]), np.mean(complex_win[:n_avg])
    else:
        amp_end = phase_end = cplx_end = amp_start = phase_start = cplx_start = 0.0

    if n_zeros > 0:
        amp_pad = np.zeros(pad_to, dtype=float)
        phase_pad = np.zeros(pad_to, dtype=float)
        complex_pad = np.zeros(pad_to, dtype=complex)
        
        if pad_mode == 'center':    
            start_idx = n_zeros // 2  
            amp_pad[:start_idx] = amp_start
            amp_pad[start_idx+n_points:] = amp_end
            phase_pad[:start_idx] = phase_start
            phase_pad[start_idx+n_points:] = phase_end
            complex_pad[:start_idx] = cplx_start
            complex_pad[start_idx+n_points:] = cplx_end
            
            amp_pad[start_idx:start_idx+n_points] = amp_win
            phase_pad[start_idx:start_idx+n_points] = phase_win
            complex_pad[start_idx:start_idx+n_points] = complex_win
            
            x_pad_pre = x_c[0] - dx * np.arange(start_idx, 0, -1)
            x_pad_post = x_c[-1] + dx * np.arange(1, (n_zeros - start_idx) + 1)
        elif pad_mode == 'zero':  # 🌟 NEW option
            amp_pad[:] = 0.0
            phase_pad[:] = 0.0
            complex_pad[:] = 0.0j
            
            amp_pad[:n_points] = amp_win
            phase_pad[:n_points] = phase_win
            complex_pad[:n_points] = complex_win
            
            x_pad_pre = []
            x_pad_post = x_c[-1] + dx * np.arange(1, n_zeros + 1)
        else:                       
            amp_pad[:] = amp_end
            phase_pad[:] = phase_end
            complex_pad[:] = cplx_end
            
            amp_pad[:n_points] = amp_win
            phase_pad[:n_points] = phase_win
            complex_pad[:n_points] = complex_win
            
            x_pad_pre = []
            x_pad_post = x_c[-1] + dx * np.arange(1, n_zeros + 1)
    else:
        amp_pad = amp_win
        phase_pad = phase_win
        complex_pad = complex_win
        x_pad_pre, x_pad_post = [], []

    pad_str = f"N={n_points}\u2192{pad_to}" if n_zeros > 0 else f"N={n_points} (no pad)"

    # 执行 FFT
    fft_amp = np.fft.rfft(amp_pad)
    q_amp = 2 * np.pi * np.fft.rfftfreq(pad_to, d=dx)
    spec_amp = np.abs(fft_amp)

    fft_phase = np.fft.rfft(phase_pad)
    spec_phase = np.abs(fft_phase)
    
    # 🌟 彻底物理剔除 q=0 直流点，防止视觉连线误导
    if len(q_amp) > 0 and q_amp[0] == 0:
        q_amp = q_amp[1:]
        spec_amp = spec_amp[1:]
        spec_phase = spec_phase[1:]
        
    if len(spec_amp) > 0 and np.max(spec_amp) > 0: spec_amp /= np.max(spec_amp)
    if len(spec_phase) > 0 and np.max(spec_phase) > 0: spec_phase /= np.max(spec_phase)
    
    fft_complex = np.fft.fftshift(np.fft.fft(complex_pad))
    q_complex_all = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(pad_to, d=dx))
    spec_complex_all = np.abs(fft_complex)
    
    mask_pos = (q_complex_all >= 0)
    q_complex = q_complex_all[mask_pos]
    spec_complex = spec_complex_all[mask_pos]
    
    # 🌟 彻底物理剔除 q=0 直流点
    if len(q_complex) > 0 and q_complex[0] == 0: 
        q_complex = q_complex[1:]
        spec_complex = spec_complex[1:]
        
    if len(spec_complex) > 0 and np.max(spec_complex) > 0: spec_complex /= np.max(spec_complex)
    
    # 获取 Peak
    def fit_peaks(qs, spec, q_guess, search_window, method):
        def find_local_max():
            extracted = []
            if q_guess is not None:
                for qg in q_guess:
                    q_mask = (qs >= qg * (1 - search_window)) & (qs <= qg * (1 + search_window))
                    if np.any(q_mask):
                        best_q = qs[q_mask][np.argmax(spec[q_mask])]
                        extracted.append(best_q)
                    else:
                        extracted.append(None)
            return {'peaks': extracted, 'fwhm': [None]*len(extracted) if extracted else []}, None
        
        if method != 'lorentzian' or not q_guess: return find_local_max()
        try:
            def multi_lorentz(x, *p):
                y = np.zeros_like(x)
                n_peaks = len(p) // 3
                for i in range(n_peaks):
                    A, x0, gamma = p[i*3:i*3+3]
                    y += A * (gamma**2 / ((x - x0)**2 + gamma**2))
                return y 
                
            p0, bounds_lower, bounds_upper = [], [], []
            for qg in q_guess:
                lower, upper = qg * (1 - search_window), qg * (1 + search_window)
                q_mask = (qs >= lower) & (qs <= upper)
                A_init = np.max(spec[q_mask]) if np.any(q_mask) else 0.5
                p0.extend([A_init, qg, qg * 0.1])      
                bounds_lower.extend([0, lower, 0.01])
                # 将 gamma 上限设为 qg (即 FWHM <= 2 * qg)，保证它至少是一个欠阻尼峰
                bounds_upper.extend([A_init * 5.0, upper, qg])
                
            popt, _ = curve_fit(multi_lorentz, qs, spec, p0=p0, bounds=(bounds_lower, bounds_upper))
            q_smooth = np.linspace(np.min(qs), np.max(qs), 400)
            components = {'q_fit': q_smooth, 'sum': multi_lorentz(q_smooth, *popt), 'offset': 0, 'peaks': []}
            
            e_peaks = []
            e_fwhm = []
            for i in range(len(q_guess)):
                A, x0, gamma = popt[i*3:i*3+3]
                e_peaks.append(x0)
                e_fwhm.append(gamma * 2.0)
                y_peak = A * (gamma**2 / ((q_smooth - x0)**2 + gamma**2))
                components['peaks'].append({'x0': x0, 'y': y_peak, 'gamma': gamma})
            return {'peaks': e_peaks, 'fwhm': e_fwhm}, components
        except Exception as e:
            return find_local_max()
            
    # q_guess 输入单位是 10⁵ cm⁻¹，内部计算用 µm⁻¹，转换系数 ×10
    q_guess_um = [qg * 10.0 for qg in q_guess] if q_guess else None
    peaks_amp,     comp_amp     = fit_peaks(q_amp,     spec_amp,     q_guess_um, search_window, peak_method)
    peaks_phase,   comp_phase   = fit_peaks(q_amp,     spec_phase,   q_guess_um, search_window, peak_method)
    peaks_complex, comp_complex = fit_peaks(q_complex, spec_complex, q_guess_um, search_window, peak_method)

    # ================= 作图 =================
    if plot:
        fig, axs = plt.subplots(1, 4, figsize=(15, 3.8), gridspec_kw={'width_ratios': [1.6, 1, 1, 1]})
        lbl_suffix = f"  {wn}" if wn else ""
        bbox_props = dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.7)
        colors = ['red', 'green', 'magenta', 'orange', 'purple', 'cyan'] 
        
        # [图1: Spatial]
        ax1 = axs[0]
        ax1_twin = ax1.twinx()
        
        if not is_folded:
            ax1.plot(x_full, amp_full, 'k-o', ms=3, lw=1.5, label=f'Amp{lbl_suffix}', zorder=3)
            ax1_twin.plot(x_full, phase_full, 'r-o', ms=3, lw=1.5, alpha=0.5, label=f'Phase', zorder=2)
            if window and window != 'boxcar':
                ax1.plot(x_c, amp_win + np.mean(amp_c), color='#1f77b4', linestyle='-', lw=2.5, alpha=0.9, label=f'Windowed Amp', zorder=4)
        else:
            # Folded mode: just plot the final folded sequence, NO original sequence.
            # Note: amp_c already contains the physical spatial mean!
            ax1.plot(x_c, amp_c, 'k-^', ms=4, lw=1.5, label=f'Folded Amp{lbl_suffix}', zorder=3)
            ax1_twin.plot(x_c, phase_c, 'r-^', ms=4, lw=1.5, alpha=0.5, label=f'Folded Phase', zorder=2)
            if window and window != 'boxcar':
                ax1.plot(x_c, amp_win + np.mean(amp_c), color='#1f77b4', linestyle='-', lw=2.5, alpha=0.9, label=f'Folded Windowed Amp', zorder=4)
                
        x_min, x_max = (np.nanmin(x_c), np.nanmax(x_c)) if is_folded else (np.nanmin(x_full), np.nanmax(x_full))
        
        # 🌟 画出 Amp 和 Phase 完整的 Padding 延长区域
        if n_zeros > 0:
            # 🌟 直接从马上要喂进 FFT 的终极数组 (amp_pad) 里抽数据来画，保证所见即真实
            if len(x_pad_post) > 0:
                post_len = len(x_pad_post)
                ax1.plot(x_pad_post, amp_pad[-post_len:] + np.mean(amp_c), 
                         color='gray', linestyle='--', lw=2, alpha=0.4, label="Pad Amp", zorder=1)
                ax1_twin.plot(x_pad_post, phase_pad[-post_len:] + np.mean(phase_c), 
                              color='r', linestyle='--', lw=2, alpha=0.2, label="Pad Phase", zorder=1)
            if len(x_pad_pre) > 0:
                pre_len = len(x_pad_pre)
                ax1.plot(x_pad_pre, amp_pad[:pre_len] + np.mean(amp_c), color='gray', linestyle='--', lw=2, alpha=0.4, zorder=1)
                ax1_twin.plot(x_pad_pre, phase_pad[:pre_len] + np.mean(phase_c), color='r', linestyle='--', lw=2, alpha=0.2, zorder=1)
                
            # 强制左侧界限
            ax1.set_xlim(np.nanmin([np.min(x_pad_pre) if len(x_pad_pre)>0 else x_min, x_min]) - 0.1, 
                         np.nanmax([np.max(x_pad_post) if len(x_pad_post)>0 else x_max, x_max]) + 0.1)
        else:
            ax1.set_xlim(x_min - 0.1, x_max + 0.1)

        # 绘制被抛弃区域的灰色背板 (只有不折叠时才需要显示全图背板隔离区)
        if not is_folded:
            valid_x_min, valid_x_max = np.min(x_c), np.max(x_c)
            if valid_x_min > x_min: ax1.axvspan(x_min, valid_x_min, color='gray', alpha=0.3, zorder=0)
            if valid_x_max < x_max: ax1.axvspan(valid_x_max, x_max, color='gray', alpha=0.3, zorder=0)

        ax1.set_xlabel(r'Distance ($\mu$m)', fontweight='bold')
        ax1.set_ylabel('Amplitude (a.u.)', color='k', fontweight='bold')
        ax1.tick_params(axis='y', labelcolor='k')
        ax1_twin.set_ylabel('Phase (a.u.)', color='r', fontweight='bold')
        ax1_twin.tick_params(axis='y', labelcolor='r')

        axs[0].set_title(f'Spatial Data ({label})  xr={xr}')
        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax1_twin.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='lower right', fontsize=9, frameon=False)

        # ── convert q to 10⁵ cm⁻¹ for display (internal q is in µm⁻¹; 1 µm⁻¹ = 10⁴ cm⁻¹ = 0.1×10⁵ cm⁻¹)
        q_amp_plt     = q_amp     / 10.0
        q_complex_plt = q_complex / 10.0
        xlabel_q = r'$q$  ($10^5$ cm$^{-1}$)'

        # [图2: Amp FFT]
        axs[1].plot(q_amp_plt, spec_amp, 'bo-', lw=1.5, ms=4)
        axs[1].set(xlabel=xlabel_q, ylabel='Norm. FFT Amp', title=f'Amp FFT ({label})\n{pad_str}')
        axs[1].set_xlabel(xlabel_q, fontweight='bold'); axs[1].set_ylabel('Norm. FFT Amp', fontweight='bold')
        axs[1].set_xlim(q_range)
        if comp_amp is not None:
            axs[1].plot(comp_amp['q_fit']/10, comp_amp['sum'], 'k--', lw=1.5, alpha=0.8, zorder=5)
            for i, pk_data in enumerate(comp_amp['peaks']):
                c = colors[i % len(colors)]
                axs[1].fill_between(comp_amp['q_fit']/10, comp_amp['offset'], comp_amp['offset'] + pk_data['y'], color=c, alpha=0.3, zorder=4)

        for i, q_max in enumerate(peaks_amp['peaks']):
            if q_max is None: continue
            c = colors[i % len(colors)]
            axs[1].axvline(q_max/10, color=c, linestyle='--', alpha=0.6)
            fwhm_str = f"\nFWHM={peaks_amp['fwhm'][i]/10:.2f}" if peaks_amp['fwhm'][i] is not None else ""
            axs[1].text(0.95, 0.95 - i*0.13, f'q={q_max/10:.2f}{fwhm_str}', color=c, fontsize=10, fontweight='bold', transform=axs[1].transAxes, ha='right', va='top', bbox=bbox_props)

        # [图3: Phase FFT]
        axs[2].plot(q_amp_plt, spec_phase, 'm^-', lw=1.5, ms=4)
        axs[2].set(xlabel=xlabel_q, ylabel='Norm. FFT Amp', title=f'Phase FFT ({label})\n{pad_str}')
        axs[2].set_xlabel(xlabel_q, fontweight='bold'); axs[2].set_ylabel('Norm. FFT Amp', fontweight='bold')
        axs[2].set_xlim(q_range)
        if comp_phase is not None:
            axs[2].plot(comp_phase['q_fit']/10, comp_phase['sum'], 'k--', lw=1.5, alpha=0.8, zorder=5)
            for i, pk_data in enumerate(comp_phase['peaks']):
                c = colors[i % len(colors)]
                axs[2].fill_between(comp_phase['q_fit']/10, comp_phase['offset'], comp_phase['offset'] + pk_data['y'], color=c, alpha=0.3, zorder=4)

        for i, q_max in enumerate(peaks_phase['peaks']):
            if q_max is None: continue
            c = colors[i % len(colors)]
            axs[2].axvline(q_max/10, color=c, linestyle='--', alpha=0.6)
            fwhm_str = f"\nFWHM={peaks_phase['fwhm'][i]/10:.2f}" if peaks_phase['fwhm'][i] is not None else ""
            axs[2].text(0.95, 0.95 - i*0.13, f'q={q_max/10:.2f}{fwhm_str}', color=c, fontsize=10, fontweight='bold', transform=axs[2].transAxes, ha='right', va='top', bbox=bbox_props)

        # [图4: Complex FFT]
        axs[3].plot(q_complex_plt, spec_complex, 'go-', lw=1.5, ms=4)
        axs[3].set(xlabel=xlabel_q, ylabel='Norm. FFT Amp', title=f'Complex FFT ({label})\n{pad_str}')
        axs[3].set_xlabel(xlabel_q, fontweight='bold'); axs[3].set_ylabel('Norm. FFT Amp', fontweight='bold')
        axs[3].set_xlim(q_range)
        if comp_complex is not None:
            axs[3].plot(comp_complex['q_fit']/10, comp_complex['sum'], 'k--', lw=1.5, alpha=0.8, zorder=5)
            for i, pk_data in enumerate(comp_complex['peaks']):
                c = colors[i % len(colors)]
                axs[3].fill_between(comp_complex['q_fit']/10, comp_complex['offset'], comp_complex['offset'] + pk_data['y'], color=c, alpha=0.3, zorder=4)

        for i, q_max in enumerate(peaks_complex['peaks']):
            if q_max is None: continue
            c = colors[i % len(colors)]
            axs[3].axvline(q_max/10, color=c, linestyle='--', alpha=0.6)
            fwhm_str = f"\nFWHM={peaks_complex['fwhm'][i]/10:.2f}" if peaks_complex['fwhm'][i] is not None else ""
            axs[3].text(0.95, 0.95 - i*0.13, f'q={q_max/10:.2f}{fwhm_str}', color=c, fontsize=10, fontweight='bold', transform=axs[3].transAxes, ha='right', va='top', bbox=bbox_props)

        for ax_i in [ax1, ax1_twin, axs[1], axs[2], axs[3]]:
            ax_i.tick_params(direction=tick_dir, which='both', top=True, right=True)
            ax_i.yaxis.set_major_locator(MaxNLocator(nbins=y_bins))

        plt.tight_layout()
        plt.show()
    
    return {
        'q_amp': q_amp, 'spec_amp': spec_amp,
        'q_phase': q_amp, 'spec_phase': spec_phase,
        'q_complex': q_complex, 'spec_complex': spec_complex,
        'peaks_amp': peaks_amp, 'peaks_phase': peaks_phase, 'peaks_complex': peaks_complex,
        'n_orig': n_points, 'n_pad': pad_to,  # expose for downstream N labeling
        # post-fold spatial data (used for spatial stack preview)
        'x_spatial': np.array(x_c),
        'y_amp_spatial': np.array(amp_c),
        'y_phase_spatial': np.array(phase_c),
    }

def plot_stacked_fft(res_dict, ch_label, wn, manual_peaks, stacked_figsize=(4, 6), 
                     stacked_xlim=(0, 5), stacked_ylim=(-0.1, 3.2), stacked_offset=0.6,
                     tick_dir='in', y_bins=5):
    fig, ax = plt.subplots(figsize=stacked_figsize)
    qa, sa = res_dict['q_amp'], res_dict['spec_amp']
    qp, sp = res_dict['q_phase'], res_dict['spec_phase']
    qc, sc = res_dict['q_complex'], res_dict['spec_complex']
    
    # convert q to 10⁵ cm⁻¹ (internal: µm⁻¹; 1 µm⁻¹ = 10⁴ cm⁻¹ = 0.1×10⁵ cm⁻¹)
    ax.plot(qa/10, sa,                    '-o', label=f'{ch_label}A', color='#4c72b0', lw=2, ms=5)
    ax.plot(qp/10, sp + stacked_offset,   '-o', label=f'{ch_label}P', color='#dd8452', lw=2, ms=5)
    ax.plot(qc/10, sc + 2*stacked_offset, '-o', label='complex',      color='#55a868', lw=2, ms=5)
    
    bbox_props = dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.8)
    colors_peak = ['red', 'blue', 'black', 'magenta', 'gray']
    
    for i, q_val in enumerate(manual_peaks):
        c = colors_peak[i % len(colors_peak)]
        # manual_peaks 单位已经是 10⁵ cm⁻¹，直接画线
        ax.axvline(q_val, color=c, linestyle='--', alpha=0.8)
        ax.text(0.3, 0.95 - i*0.08, f'q$_{{{i+1}}}$ = {q_val:.2f}', 
                color=c, fontsize=12, fontweight='bold',
                transform=ax.transAxes, ha='left', va='top', bbox=bbox_props)
    
    ax.set_xlim(stacked_xlim)
    if stacked_ylim:
        ax.set_ylim(stacked_ylim)
        
    ax.set_xlabel(r'$q$  ($10^5$ cm$^{-1}$)', fontsize=12)
    ax.set_ylabel('FFT Amplitude (a.u.)', fontsize=12)
    ax.set_title(f'Stacked FFT for {ch_label} ({wn})')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.tick_params(direction=tick_dir, which='both')
    ax.yaxis.set_major_locator(MaxNLocator(nbins=y_bins))
    
    ax.legend(loc='upper right', frameon=False, fontsize=11)
    plt.tight_layout()
    plt.show()

def analyze_channel_ffts(amplp, phaselp, wn, 
                         xr_O1=(0.3, 3.5), xr_O2=(0.28, 3.5),
                         mirror_at_O1=None, mirror_at_O2=None,
                         q_range=(0, 5),                             # 10⁵ cm⁻¹, all FFT panels
                         window='boxcar', pad_factor=2.0, pad_mode='end',
                         q_guess=[5, 20], search_window=0.3, peak_method='find_peak',
                         tick_dir='in', y_bins=5,
                         stacked_figsize=(4, 6), stacked_ylim=(-0.1, 3.2), stacked_offset=0.6,
                         manual_peaks_O1=[5, 20], manual_peaks_O2=[5, 20],
                         plot=True, plot_stacked=True):
    """
    Main wrapper to run the Stacked FFT analysis on both O1 and O2 channels.
    Takes amplp and phaselp Line Profiles DataFrames and a wavenumber string (e.g. '340.00cm-1').
    Returns a dictionary with comprehensive FFT data and extracted peaks.
    """
    import re
    
    # 🌟 智能寻列：通过用户给的 wn (例如 '340.00cm-1' 或 '340') 在 df 里模糊匹配出真实的列名
    match_target = re.search(r'[\d\.]+', str(wn))
    target_val = float(match_target.group()) if match_target else 0.0
    
    real_wns = []
    # 搜集数据表中的所有真实频率
    for col in amplp.columns:
        if col.endswith('_O1A') or col.endswith('_O2A'):
            match = re.match(r'^([\d\.]+)cm\-1_', col)
            if match:
                real_val = float(match.group(1))
                real_str = col.split('_')[0]  # 形如 '340.83cm-1'
                real_wns.append((real_val, real_str))
                
    actual_wn_str = str(wn) # default fallback
    if real_wns:
        # 寻找差值最小的最接近频率
        closest = min(real_wns, key=lambda x: abs(x[0] - target_val))
        actual_wn_str = closest[1]
        print(f"🔍 智能匹配波数: 你的索引 '{wn}' -> 成功匹配到真实切片 '{actual_wn_str}'")
        
    col_O1A, col_O1P = f'{actual_wn_str}_O1A', f'{actual_wn_str}_O1P'
    col_O2A, col_O2P = f'{actual_wn_str}_O2A', f'{actual_wn_str}_O2P'
    
    x_dist = amplp['distance_um']
    
    # 🌟 返回的字典会同时保存你的 'index'(wn_index) 和 '真实值'(actual_wn)
    results = {'wn_index': wn, 'actual_wn': actual_wn_str, 'channels': {}}
    
    if col_O1A in amplp.columns and col_O1P in phaselp.columns:
        if plot:
            print(f">>> Drawing Channel O1 for {actual_wn_str} ...")
        # 作图时传入真实波数 actual_wn_str，确保图表 Title 和 Label 精准备注
        res_O1 = plot_channel_fft(x_dist, amplp[col_O1A], phaselp[col_O1P], 
                                  'O1', wn=actual_wn_str, xr=xr_O1, mirror_at=mirror_at_O1,
                                  q_range=q_range,
                                  window=window, pad_factor=pad_factor, pad_mode=pad_mode, 
                                  q_guess=q_guess, search_window=search_window,
                                  peak_method=peak_method, tick_dir=tick_dir, y_bins=y_bins, plot=plot)
        
        if plot:
            print("\n=== ✨ Extracted Peaks for O1 ===")
            print("Amp FFT:     ", res_O1['peaks_amp'])
            print("Phase FFT:   ", res_O1['peaks_phase'])
            print("Complex FFT: ", res_O1['peaks_complex'])
            print("-" * 40)
            
            if plot_stacked:
                print(f"\n>>> Drawing Stacked Plot for O1 ({actual_wn_str})...")
                plot_stacked_fft(res_O1, 'O1', actual_wn_str, manual_peaks_O1, stacked_figsize=stacked_figsize,
                                 stacked_xlim=q_range, stacked_ylim=stacked_ylim, stacked_offset=stacked_offset,
                                 tick_dir=tick_dir, y_bins=y_bins)
                         
        results['channels']['O1'] = res_O1
    else:
        print(f"⚠️ Warning: O1 columns ({col_O1A}, {col_O1P}) not found in data.")

    if col_O2A in amplp.columns and col_O2P in phaselp.columns:
        if plot:
            print(f">>> Drawing Channel O2 for {actual_wn_str} ...")
        res_O2 = plot_channel_fft(x_dist, amplp[col_O2A], phaselp[col_O2P], 
                                  'O2', wn=actual_wn_str, xr=xr_O2, mirror_at=mirror_at_O2,
                                   q_range=q_range,
                                  window=window, pad_factor=pad_factor, pad_mode=pad_mode, 
                                  q_guess=q_guess, search_window=search_window,
                                  peak_method=peak_method, tick_dir=tick_dir, y_bins=y_bins, plot=plot)
        
        if plot:
            print("\n=== ✨ Extracted Peaks for O2 ===")
            print("Amp FFT:     ", res_O2['peaks_amp'])
            print("Phase FFT:   ", res_O2['peaks_phase'])
            print("Complex FFT: ", res_O2['peaks_complex'])
            print("=" * 40)
            
            if plot_stacked:
                print(f"\n>>> Drawing Stacked Plot for O2 ({actual_wn_str})...")
                plot_stacked_fft(res_O2, 'O2', actual_wn_str, manual_peaks_O2, stacked_figsize=stacked_figsize,
                                 stacked_xlim=q_range, stacked_ylim=stacked_ylim, stacked_offset=stacked_offset,
                                 tick_dir=tick_dir, y_bins=y_bins)
                         
        results['channels']['O2'] = res_O2
    else:
        print(f"⚠️ Warning: O2 columns ({col_O2A}, {col_O2P}) not found in data.")
        
    return results

def batch_analyze_and_plot_stacked_ffts(amplp, phaselp,
                                        xr_O1=(0.27, 3.5), xr_O2=(0.27, 3.5),
                                        mirror_at_O1=None, mirror_at_O2=None,
                                        window='boxcar', pad_factor=3, pad_mode='end',
                                        q_guess=[8, 19], search_window=0.2, peak_method='find_peak',
                                        q_range=(0, 25), stacked_figsize=(11, 7.5), stacked_offset=0.25,
                                        save_dir='data/NearField/lineprofile/Sample11'):
    """
    Batches through all available wavenumbers in the dataframes, performs FFT analysis exactly like
    `analyze_channel_ffts`, but generates a global 3-panel stacked plot (Amplitude, Phase, Complex)
    for both O1 and O2 (just like typical waterfall maps), and automates saving the results to a CSV
    whose name reflects the configured padding and mirror settings.
    """
    import re
    import os
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    import matplotlib.cm as cm
    import sys
    import numpy as np
    
    os.makedirs(save_dir, exist_ok=True)

    # 1. Gather all actual wavenumbers from the dataframe
    real_wns = []
    for col in amplp.columns:
        if col.endswith('_O1A') or col.endswith('_O2A'):
            match = re.match(r'^([\d\.]+)cm\-1_', col)
            if match:
                real_val = float(match.group(1))
                real_str = col.split('_')[0]
                if (real_val, real_str) not in real_wns:
                    real_wns.append((real_val, real_str))
                    
    # Sort them by frequency
    real_wns.sort(key=lambda x: x[0])
    
    if not real_wns:
        print("❌ No matching columns found in amplp/phaselp!")
        return {}

    print(f"🚀 Batch processing {len(real_wns)} wavenumbers...")
    
    all_fft_results = {}
    
    # 2. Run FFT for each (silently!)
    for idx, (wn_val, wn_str) in enumerate(real_wns):
        # --- 🌟 智能动态 `q_guess` 逻辑 ---
        # 支持三种模式：
        #   - 全局模式：q_guess=[8,19]，长度不等于切片数 → 所有切片都用这个
        #   - 单峰追踪：q_guess=[4,5,6,None,8,...]，长度等于切片数，标量→单峰，None→跳过
        #   - 多峰追踪：q_guess=[[8,19],[8.5,19.5],...]，每个元素是列表
        current_q_guess = q_guess
        if isinstance(q_guess, (list, tuple)) and len(q_guess) == len(real_wns):
            elem = q_guess[idx]
            if elem is None:
                # None 表示这个切片没有可信的峰，跳过标注
                current_q_guess = None
            elif hasattr(elem, '__iter__') and not isinstance(elem, str):
                # 嵌套列表：[[8,19], [8.5,19.5], ...]
                current_q_guess = [x for x in elem if x is not None]
            else:
                # 标量：强转为单元素列表 → 只找一个峰
                current_q_guess = [elem]
                
        res = analyze_channel_ffts(
            amplp, phaselp, wn=wn_str, 
            xr_O1=xr_O1, xr_O2=xr_O2,
            mirror_at_O1=mirror_at_O1, mirror_at_O2=mirror_at_O2,
            q_range_amp=q_range, q_range_complex=q_range,
            window=window, pad_factor=pad_factor, pad_mode=pad_mode,
            q_guess=current_q_guess, search_window=search_window, peak_method=peak_method,
            plot=False  # Important!
        )
        all_fft_results[wn_str] = res
        peak_info = f"q≈{current_q_guess}" if current_q_guess else "no peak"
        print(f"  ✓ {wn_str}  ({peak_info})      ", end='\r')

        
    print("\n✅ All FFT processing compelted. Generating global stacked maps...")
    
    # 3. Generating Plot for each Channel
    def draw_spatial_stack(ch_label):
        """Stacked real-space (post-fold) linecuts for all wavenumbers."""
        valid_res = [(wn_val, wn_str, all_fft_results[wn_str]) for wn_val, wn_str in real_wns
                     if ch_label in all_fft_results.get(wn_str, {}).get('channels', {})]
        if not valid_res:
            return
        _r0 = valid_res[0][2]['channels'][ch_label]
        _n_orig = _r0.get('n_orig', 0)
        _n_pad  = _r0.get('n_pad', _n_orig)
        _npad_tag = f"N={_n_orig}\u2192{_n_pad}(\u00d7{pad_factor})" if _n_pad > _n_orig else f"N={_n_orig}"
        _mir = mirror_at_O1 if ch_label == 'O1' else mirror_at_O2

        fig_sp, ax_sp = plt.subplots(figsize=(5.5, stacked_figsize[1]), constrained_layout=True)
        fig_sp.suptitle(f"{ch_label} Real-space stacked linecuts  [{_npad_tag}]  mirror={_mir}", fontsize=10)
        cmap_sp = cm.get_cmap('Spectral_r')
        wn_min_sp, wn_max_sp = valid_res[0][0], valid_res[-1][0]
        _sp_offset = stacked_offset * 0.35

        for idx, (wn_val, wn_str, master_res) in enumerate(valid_res):
            rd = master_res['channels'][ch_label]
            x_s = rd.get('x_spatial')
            y_s = rd.get('y_amp_spatial')
            if x_s is None or y_s is None or len(x_s) == 0:
                continue
            c = cmap_sp((wn_val - wn_min_sp) / max(wn_max_sp - wn_min_sp, 1e-5))
            y_base = idx * _sp_offset
            y_norm = y_s - np.mean(y_s)
            y_rng = np.max(np.abs(y_norm)) if np.max(np.abs(y_norm)) > 0 else 1.0
            ax_sp.plot(x_s, y_norm / y_rng * (_sp_offset * 0.45) + y_base, '-', lw=1.2, color=c)
            ax_sp.text(x_s[-1] + 0.05, y_base, f"{wn_val:.0f}", color=c, fontsize=7, va='center')

        ax_sp.set_xlabel(r'Distance from edge ($\mu$m)')
        ax_sp.set_yticks([])
        ax_sp.spines['top'].set_visible(False)
        ax_sp.spines['right'].set_visible(False)
        ax_sp.spines['left'].set_visible(False)
        plt.show()

    def draw_global_stack(ch_label):
        """Stacked FFT waterfall with per-channel N label."""
        valid_res = [(wn_val, wn_str, all_fft_results[wn_str]) for wn_val, wn_str in real_wns
                     if ch_label in all_fft_results.get(wn_str, {}).get('channels', {})]
        if not valid_res:
            print(f"\u26a0\ufe0f Skipping {ch_label} stacked plot (No data)")
            return
        _first = valid_res[0][2]['channels'][ch_label]
        _n_orig = _first.get('n_orig', 0)
        _n_pad  = _first.get('n_pad', _n_orig)
        _npad_tag = f"N={_n_orig}\u2192{_n_pad}(\u00d7{pad_factor})" if _n_pad > _n_orig else f"N={_n_orig}"
        _mir = mirror_at_O1 if ch_label == 'O1' else mirror_at_O2
        title_params = f"mirror={_mir}  [{_npad_tag}]  {pad_mode}-pad"

        fig, axs = plt.subplots(1, 3, figsize=stacked_figsize, sharey=False)
        fig.suptitle(f"{ch_label} stacked 1-D FFT  |  {title_params}", fontsize=11, y=0.98)

        ax_amp, ax_phase, ax_complex = axs
        ax_amp.set_title(f"{ch_label} Amplitude", pad=15)
        ax_phase.set_title(f"{ch_label} Phase", pad=15)
        ax_complex.set_title(f"{ch_label} |A+iP|", pad=15)

        cmap = cm.get_cmap('Spectral_r')
        wn_min, wn_max = valid_res[0][0], valid_res[-1][0]
        
        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=wn_min, vmax=wn_max))
        sm.set_array([])
        cbar_ax = fig.add_axes([0.92, 0.35, 0.015, 0.3])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Frequency (cm$^{-1}$)')
        cbar.ax.yaxis.set_ticks_position('right')
        cbar.ax.yaxis.set_label_position('right')

        for idx, (wn_val, wn_str, master_res) in enumerate(valid_res):
            res_dict = master_res['channels'][ch_label]
            c = cmap((wn_val - wn_min) / max(wn_max - wn_min, 1e-5))
            y_base = idx * stacked_offset
            
            qa, sa = res_dict['q_amp'], res_dict['spec_amp']
            qp, sp = res_dict['q_phase'], res_dict['spec_phase']
            qc, sc = res_dict['q_complex'], res_dict['spec_complex']
            
            # Plot traces
            ax_amp.plot(qa, sa + y_base, '.-', ms=4, color=c, lw=1.2)
            ax_phase.plot(qp, sp + y_base, '.-', ms=4, color=c, lw=1.2)
            ax_complex.plot(qc, sc + y_base, '.-', ms=4, color=c, lw=1.2)
            
            # Mark peaks (Amplitude)
            for i, p_val in enumerate(res_dict['peaks_amp']['peaks']):
                if p_val is not None:
                    idx_q = np.argmin(np.abs(qa - p_val))
                    ax_amp.plot(p_val, sa[idx_q] + y_base, 'o', mec='#4c72b0', mfc='none', mew=1.5, ms=7)
                    txt_va = 'bottom' if idx % 2 == 0 else 'top'
                    txt_dy = +0.05 if idx % 2 == 0 else -0.05
                    ax_amp.text(p_val, sa[idx_q] + y_base + txt_dy, f"{p_val:.1f}",
                                color='#4c72b0', fontsize=8, ha='center', va=txt_va)
                                    
            # Peak marking for Phase
            for i, p_val in enumerate(res_dict['peaks_amp']['peaks']):
                if p_val is not None:
                    idx_q = np.argmin(np.abs(qp - p_val))
                    ax_phase.plot(p_val, sp[idx_q] + y_base, 'o', mec='#dd8452', mfc='none', mew=1.5, ms=7)
                    txt_va = 'bottom' if idx % 2 == 0 else 'top'
                    txt_dy = +0.05 if idx % 2 == 0 else -0.05
                    ax_phase.text(p_val, sp[idx_q] + y_base + txt_dy, f"{p_val:.1f}",
                                  color='#dd8452', fontsize=8, ha='center', va=txt_va)

            # Mark peaks (Complex)
            for i, p_val in enumerate(res_dict['peaks_complex']['peaks']):
                if p_val is not None:
                    idx_q = np.argmin(np.abs(qc - p_val))
                    ax_complex.plot(p_val, sc[idx_q] + y_base, 'o', mec='#55a868', mfc='none', mew=1.5, ms=7)
                    txt_va = 'bottom' if idx % 2 == 0 else 'top'
                    txt_dy = +0.05 if idx % 2 == 0 else -0.05
                    ax_complex.text(p_val, sc[idx_q] + y_base + txt_dy, f"{p_val:.1f}",
                                    color='#55a868', fontsize=8, ha='center', va=txt_va)

            # Annotate frequency tick slightly off to the right
            for ax_i in axs:
                ax_i.text(q_range[1] + 1, y_base + 0.1, f"{wn_val:.0f}", color=c, fontsize=7, alpha=0.8, va='center')

        # Global Formatting
        for ax in axs:
            ax.set_xlim(q_range)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.set_yticks([])
            ax.set_xlabel(r'q ($\mu$m$^{-1}$)')
        
        plt.subplots_adjust(wspace=0.1)
        plt.show()

    # ── Draw: spatial stack + FFT waterfall for each channel ──
    draw_spatial_stack('O1')
    draw_global_stack('O1')
    draw_spatial_stack('O2')
    draw_global_stack('O2')

    # 4. Save to CSV
    import pandas as pd
    rows = []
    for wn_index, data in all_fft_results.items():
        row_dict = {
            'wn_index': wn_index,
            'actual_wn': data['actual_wn']
        }
        for ch in ['O1', 'O2']:
            if ch in data['channels']:
                ch_data = data['channels'][ch]
                for i, (peak, fwhm) in enumerate(zip(ch_data['peaks_amp']['peaks'], ch_data['peaks_amp']['fwhm'])):
                    row_dict[f'{ch}_Amp_peak_{i+1}'] = peak
                    row_dict[f'{ch}_Amp_fwhm_{i+1}'] = fwhm
                for i, (peak, fwhm) in enumerate(zip(ch_data['peaks_complex']['peaks'], ch_data['peaks_complex']['fwhm'])):
                    row_dict[f'{ch}_Cplx_peak_{i+1}'] = peak
                    row_dict[f'{ch}_Cplx_fwhm_{i+1}'] = fwhm
        rows.append(row_dict)

    df_fft_export = pd.DataFrame(rows)
    # Dynamic filename avoiding invalid chars
    def safe_str(val): return str(val).replace('.', 'p')
    csv_name = f"Extracted_Peaks_pad{safe_str(pad_factor)}_{pad_mode}_mO1_{safe_str(mirror_at_O1)}_mO2_{safe_str(mirror_at_O2)}.csv"
    csv_path = os.path.join(save_dir, csv_name)
    df_fft_export.to_csv(csv_path, index=False)
    
    print(f"\n✅ All results successfully saved to: \n📁 {csv_path}\n")
    return all_fft_results

# ==========================================
# 🌀 Complex Hankel Transform (CHT) Fitting
# ==========================================
def complex_hankel_transform(x, signal, L=1.0, k_array=None):
    r"""
    Performs the Complex Hankel Transform (Woessner et al. Nat Mater 2015).
    T(k) = 0.5 * int_0^L x [H_0^{(1)}(kx)]^* \delta\xi_{opt}(x) u(x) dx
    """
    import numpy as np
    import scipy.special as sp
    
    if k_array is None:
        k_array = np.linspace(0.1, 40, 400)
        
    u = np.zeros_like(x)
    mask = x <= L
    u[mask] = 1 - np.sin(np.pi / 2 * x[mask] / L)**2
    
    dx = np.gradient(x) if len(x) > 1 else np.ones_like(x)
    T_k = np.zeros(len(k_array), dtype=complex)
    
    for i, k in enumerate(k_array):
        x_safe = np.maximum(x, 1e-5)
        H0_conj = np.conj(sp.hankel1(0, k * x_safe))
        integrand = x * H0_conj * signal * u
        T_k[i] = 0.5 * np.sum(integrand * dx)
        
    return k_array, T_k

def fit_cht_peaks(x, signal, L=1.0, k_fit_range=(0.5, 30), p0=None, bounds=None, k_plot_range=(0.1, 40), num_peaks=1, model_type='independent'):
    """
    Fits the CHT peaks to extract q_p.
    Supports model_type='independent' or 'linked_paper'.
    """
    import numpy as np
    from scipy.optimize import least_squares
    import scipy.special as sp
    
    k_array_full = np.linspace(k_plot_range[0], k_plot_range[1], 400)
    _, T_data_full = complex_hankel_transform(x, signal, L=L, k_array=k_array_full)
    
    mask = (k_array_full >= k_fit_range[0]) & (k_array_full <= k_fit_range[1])
    k_array_fit = k_array_full[mask]
    T_data_fit = T_data_full[mask]
    
    def model_signal(x_val, theta, B_offset=0):
        sig = np.zeros_like(x_val, dtype=complex)
        x_shifted = np.maximum(x_val, 1e-5)
        
        if model_type == 'independent':
            for i in range(num_peaks):
                A = theta[4*i]
                q_re = theta[4*i+1]
                q_im = theta[4*i+2]
                phase = theta[4*i+3]
                q_p = q_re + 1j * q_im
                sig += A * np.exp(1j * phase) * sp.hankel1(0, 2 * q_p * x_shifted)
        elif model_type == 'linked_paper':
            A = theta[0]
            q_re = theta[1]
            q_im = theta[2]
            phase_A = theta[3]
            B_amp = theta[4]
            phase_B = theta[5]
            q_p = q_re + 1j * q_im
            
            term1 = A * np.exp(1j * phase_A) * sp.hankel1(0, 2 * q_p * x_shifted)
            term2 = B_amp * np.exp(1j * phase_B) * np.exp(1j * q_p * x_shifted) / np.sqrt(x_shifted)
            sig += term1 + term2
            
        return B_offset + sig

    def resid(theta):
        mod_sig = model_signal(x, theta)
        _, T_mod_fit = complex_hankel_transform(x, mod_sig, L=L, k_array=k_array_fit)
        diff = np.sqrt(k_array_fit) * (T_data_fit - T_mod_fit)
        return np.concatenate([np.real(diff), np.imag(diff)])
        
    if p0 is None:
        if model_type == 'independent':
            p0 = []
            for i in range(num_peaks):
                p0.extend([np.nanmax(np.abs(signal))/(i+1), 10.0 * (i+1), 0.5, 0.0])
        elif model_type == 'linked_paper':
            p0 = [np.nanmax(np.abs(signal))/2.0, 10.0, 0.5, 0.0, np.nanmax(np.abs(signal))/2.0, 0.0]
            
    if bounds is None:
        if model_type == 'independent':
            lower = [0, 0.1, 0.01, -np.pi] * num_peaks
            upper = [np.inf, 100.0, 10.0, np.pi] * num_peaks
            bounds = (lower, upper)
        elif model_type == 'linked_paper':
            lower = [0, 0.1, 0.01, -np.pi, 0.0, -np.pi]
            upper = [np.inf, 100.0, 10.0, np.pi, np.inf, np.pi]
            bounds = (lower, upper)
        
    res = least_squares(resid, p0, bounds=bounds, loss='soft_l1')
    
    if model_type == 'independent':
        fit_params = []
        for i in range(num_peaks):
            fit_params.append({
                'A': res.x[4*i],
                'q_re': res.x[4*i+1],
                'q_im': res.x[4*i+2],
                'phase': res.x[4*i+3]
            })
    elif model_type == 'linked_paper':
        fit_params = [{
            'A': res.x[0],
            'q_re': res.x[1],
            'q_im': res.x[2],
            'phase_A': res.x[3],
            'B_amp': res.x[4],
            'phase_B': res.x[5]
        }]
        
    mod_sig_fit = model_signal(x, res.x)
    _, T_mod_full = complex_hankel_transform(x, mod_sig_fit, L=L, k_array=k_array_full)

    return fit_params, k_array_full, T_data_full, T_mod_full, mod_sig_fit


def load_aligned_wn_signal(wn, align_shift_nm, data_dir='data/graphene_3x1', L_cutoff_fit=1.5):
    """
    Load one wavenumber's averaged line-profile CSV, shift its edge to x=0 using
    align_shift_nm, and background-subtract O3A with the same savgol window used
    throughout fitting_pipeline.ipynb. Returns the raw aligned dataframe plus the
    masked/background-subtracted (x_f, sig_f) used by the CHT and real-space fits.
    """
    import glob, re
    file_paths = glob.glob(f'{data_dir}/*_{wn.replace("cm-1","")}*AVG_lp1.csv')
    if not file_paths:
        # fall back to scanning all files and matching the leading wavenumber
        file_paths = [p for p in glob.glob(f'{data_dir}/*_AVG_lp1.csv')
                      if re.search(rf'(?<!\d){re.escape(wn.replace("cm-1",""))}cm-1', p)]
    if not file_paths:
        raise FileNotFoundError(f"No data file found for {wn} in {data_dir}")

    df_target = pd.read_csv(file_paths[0])
    df_target['distance_um'] = (df_target['distance_nm'] - align_shift_nm) / 1000.0

    x_mat = df_target['distance_um'].values
    y_mat = df_target['O3A'].values
    window_len = min(41, len(y_mat) if len(y_mat) % 2 != 0 else len(y_mat) - 1)
    y_bg = savgol_filter(y_mat, window_length=window_len, polyorder=2)
    y_osc = y_mat - y_bg

    mask_fit = (x_mat >= 0) & (x_mat <= L_cutoff_fit)
    x_f = x_mat[mask_fit]
    sig_f = y_osc[mask_fit]

    return dict(df_target=df_target, x_f=x_f, sig_f=sig_f)


def fit_and_plot_cht(x_f, sig_f, wn, x_start_cht=0.22, L_cutoff_cht=1.2,
                      k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=2.0,
                      model_type='linked_paper', num_peaks=1, save_path=None):
    """
    Single-source-of-truth CHT fit + 3-panel plot, parametrized version of the
    'Fit 1: Complex Hankel Transform' cell that used to be copy-pasted once per
    wavenumber in fitting_pipeline.ipynb. Same math, same plot layout; only the
    tunable parameters (k_fit_range_cm, k_linked_guess_cm, x_start_cht, ...) are
    now arguments instead of being retyped in every cell.

    Returns a dict with the fit results (lambda_p, damping, k-space relative RMSE)
    and the matplotlib figure. If save_path is given, the figure is saved there.
    """
    k_fit_range = (k_fit_range_cm[0] * 10.0, k_fit_range_cm[1] * 10.0)

    mask_cht = x_f >= x_start_cht
    x_f_cht = x_f[mask_cht]
    sig_f_cht = sig_f[mask_cht]

    amp_max = np.nanmax(np.abs(sig_f_cht))
    if model_type == 'linked_paper':
        q_re_guess = k_linked_guess_cm * 10.0
        p0_guess = [amp_max / 2.0, q_re_guess, 0.5, 0.0, amp_max / 2.0, 0.0]
    else:
        raise NotImplementedError("fit_and_plot_cht currently only wraps model_type='linked_paper'; "
                                   "use fit_cht_peaks directly for 'independent'.")

    fit_params_cht, k_arr, T_data, T_mod, mod_sig_cplx = fit_cht_peaks(
        x_f_cht, sig_f_cht, L=L_cutoff_cht, k_fit_range=k_fit_range, p0=p0_guess,
        k_plot_range=(0.1, 100), num_peaks=num_peaks, model_type=model_type
    )

    params = fit_params_cht[0]
    q_re_fit, q_im_fit = params['q_re'], params['q_im']
    q_p = q_re_fit + 1j * q_im_fit
    lam_cht_nm = (2 * np.pi / q_re_fit) * 1000
    damping_cht = q_re_fit / q_im_fit if q_im_fit > 1e-9 else np.inf

    fit_mask = (k_arr >= k_fit_range[0]) & (k_arr <= k_fit_range[1])
    rmse_k = np.sqrt(np.mean(np.abs(T_data[fit_mask] - T_mod[fit_mask]) ** 2))
    rel_rmse_k = rmse_k / np.mean(np.abs(T_data[fit_mask]))

    import scipy.special as sp
    x_safe = np.maximum(x_f_cht, 1e-5)
    env_tip = params['A'] * np.abs(sp.hankel1(0, 2 * q_p * x_safe))
    env_edge = params['B_amp'] * np.abs(np.exp(1j * q_p * x_safe) / np.sqrt(x_safe))
    envelope_list = [env_tip, env_edge]
    envelope = envelope_list[0]
    mod_sig_fit = np.real(mod_sig_cplx)

    # Dense grid purely for smooth plotted curves (envelopes + total fit). The
    # actual fit/RMSE above is unchanged, still evaluated on the real data grid.
    x_dense = np.linspace(x_start_cht, L_cutoff_cht, 600)
    x_dense_safe = np.maximum(x_dense, 1e-5)
    x_nm_dense = x_dense * 1000
    phase_A = params.get('phase_A', params.get('phase', 0.0))
    mod_sig_dense = np.real(
        params['A'] * np.exp(1j * phase_A) * sp.hankel1(0, 2 * q_p * x_dense_safe)
        + params['B_amp'] * np.exp(1j * params['phase_B']) * np.exp(1j * q_p * x_dense_safe) / np.sqrt(x_dense_safe)
    )
    env_tip_dense = params['A'] * np.abs(sp.hankel1(0, 2 * q_p * x_dense_safe))
    env_edge_dense = params['B_amp'] * np.abs(np.exp(1j * q_p * x_dense_safe) / np.sqrt(x_dense_safe))
    envelope_list_dense = [env_tip_dense, env_edge_dense]
    envelope_dense = envelope_list_dense[0]

    c_data, c_fit = '#555555', '#b2182b'
    c_env_a, c_env_b_blue, c_env_b_red = '#fddbc7', '#92c5de', '#f4a582'
    colors_peaks, colors_text = ['#d6604d', 'darkorange'], ['#b2182b', '#b05b00']
    labels_linked = ['Tip-Launched (2q)', 'Edge-Launched (q)']

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), dpi=100)

    ax = axes[0]
    x_nm = x_f_cht * 1000
    ax.fill_between(x_nm_dense, envelope_dense, -envelope_dense, color=c_env_a, alpha=0.6)
    for i, env_i in enumerate(envelope_list_dense):
        ax.plot(x_nm_dense, env_i, color=colors_peaks[i], linestyle=':', lw=1.5, label=labels_linked[i])
        ax.plot(x_nm_dense, -env_i, color=colors_peaks[i], linestyle=':', lw=1.5)
    ax.plot(x_f * 1000, sig_f, marker='x', markersize=7, markeredgewidth=1.5, linestyle='None',
            color='lightgray', label='Discarded Data')
    ax.plot(x_nm, sig_f_cht, marker='x', markersize=7, markeredgewidth=1.5, linestyle='None',
            color=c_data, label='Fitted Data')
    ax.plot(x_nm_dense, mod_sig_dense, color=c_fit, lw=2, label='Total Fit')
    ax.set_xlabel('Distance from edge (nm)', fontweight='bold')
    ax.set_ylabel(r'Re $\xi_{\mathbf{opt}}$ (a.u.)', fontweight='bold')
    ax.set_xlim(-20, L_cutoff_cht * 1000 + 20)
    ax.set_yticks([])
    ax.text(0.95, 0.95, rf"Single Mode $\lambda_p = {lam_cht_nm:.1f}$ nm",
            transform=ax.transAxes, fontsize=12, va='top', ha='right', color=colors_text[0])
    ax.legend(loc='lower right', fontsize=11, frameon=False)
    ax.tick_params(direction='in', top=True, right=True)

    ax = axes[1]
    sqrt_x = np.sqrt(x_safe)
    sqrt_x_dense = np.sqrt(x_dense_safe)
    q_im_ideal = q_re_fit / 70.0
    q_p_ideal = q_re_fit + 1j * q_im_ideal
    envelope_ideal_dense = params['A'] * np.abs(sp.hankel1(0, 2 * q_p_ideal * x_dense_safe))
    envelope_sqrt_ideal_dense = envelope_ideal_dense * sqrt_x_dense
    ax.fill_between(x_nm_dense, envelope_sqrt_ideal_dense, -envelope_sqrt_ideal_dense, color=c_env_b_blue, alpha=0.8)
    ax.fill_between(x_nm_dense, envelope_dense * sqrt_x_dense, -envelope_dense * sqrt_x_dense, color=c_env_b_red, alpha=0.8)
    ax.plot(x_nm_dense, envelope_sqrt_ideal_dense, color='#053061', linestyle=':', lw=1.5)
    ax.plot(x_nm_dense, -envelope_sqrt_ideal_dense, color='#053061', linestyle=':', lw=1.5)
    for i, env_i in enumerate(envelope_list_dense):
        env_sqrt_i = env_i * sqrt_x_dense
        ax.plot(x_nm_dense, env_sqrt_i, color=colors_peaks[i], linestyle='--', lw=1.5)
        ax.plot(x_nm_dense, -env_sqrt_i, color=colors_peaks[i], linestyle='--', lw=1.5)
    ax.plot(x_nm, sig_f_cht * sqrt_x, marker='x', markersize=7, markeredgewidth=1.5, linestyle='None', color=c_data)
    ax.plot(x_nm_dense, mod_sig_dense * sqrt_x_dense, color=c_fit, lw=2)
    ax.set_xlabel('Distance from edge (nm)', fontweight='bold')
    ax.set_ylabel(r'Re $\xi_{\mathbf{opt}} \times \sqrt{\mathbf{x}}$ (a.u.)', fontweight='bold')
    ax.set_xlim(-20, L_cutoff_cht * 1000 + 20)
    ax.set_yticks([])
    for i, env_i in enumerate(envelope_list_dense):
        lbl = "Tip-launched" if i == 0 else "Edge-launched"
        env_sqrt_i = env_i * sqrt_x_dense
        max_env = np.max(env_sqrt_i)
        y_pos = env_sqrt_i[-1]
        y_text_offset = max_env * (0.8 + 0.5 * i) if i % 2 == 0 else -max_env * (0.8 + 0.5 * i)
        ax.annotate(rf"{lbl} $\gamma_p^{{-1}} = {damping_cht:.1f}$",
                    xy=(x_nm_dense[-1] * 0.7, y_pos), xytext=(x_nm_dense[-1] * 0.15, y_pos + y_text_offset),
                    color=colors_text[i], fontsize=12, fontweight='bold',
                    arrowprops=dict(arrowstyle="->", color=colors_text[i], lw=1.5, linestyle='--'))
    ax.annotate(r"Ideal $\gamma_p^{-1} = 70$",
                xy=(x_nm_dense[-1] * 0.7, -envelope_sqrt_ideal_dense[-1]),
                xytext=(x_nm_dense[-1] * 0.15, -envelope_sqrt_ideal_dense[-1] - 0.5 * np.max(envelope_dense * sqrt_x_dense)),
                color='#053061', fontsize=12, fontweight='bold',
                arrowprops=dict(arrowstyle="->", color='#053061', lw=1.5, linestyle=':'))
    ax.tick_params(direction='in', top=True, right=True)

    ax = axes[2]
    k_arr_cm = k_arr / 10.0
    ax.plot(k_arr_cm, np.abs(T_data), 'o', color=c_data, markersize=4, label='Data $|T(k)|$')
    ax.plot(k_arr_cm, np.abs(T_mod), '-', color=c_fit, lw=2, label='Total Fit $|T(k)|$')
    ax.axvspan(k_fit_range_cm[0], k_fit_range_cm[1], color='gray', alpha=0.2, label='Fit Range')
    ax.axvline(q_re_fit / 10.0, color='purple', linestyle='--', lw=1, label='q (Edge-launched)')
    ax.axvline(2 * q_re_fit / 10.0, color='blue', linestyle='--', lw=1, label='2q (Tip-launched)')
    ax.text(0.95, 0.4, rf"Single $q_p = {q_re_fit/10.0:.2f} + i{q_im_fit/10.0:.2f}$",
            transform=ax.transAxes, fontsize=11, va='center', ha='right', color='black')
    ax.set_xlabel(r'Momentum $k$ ($10^5$ cm$^{-1}$)', fontweight='bold')
    ax.set_ylabel(r'$|T(k)|$ (a.u.)', fontweight='bold')
    ax.set_xlim(0, 10)
    ax.legend(loc='upper right', fontsize=10, frameon=False)
    ax.tick_params(direction='in', top=True, right=True)

    fig.suptitle(f"CHT Fit for {wn}", fontsize=14, fontweight='bold')
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches='tight')

    results = dict(
        wn=wn, q_re=q_re_fit, q_im=q_im_fit, lambda_p_nm=lam_cht_nm,
        damping=damping_cht, rel_rmse_k=rel_rmse_k,
        k_fit_range_cm=k_fit_range_cm, k_linked_guess_cm=k_linked_guess_cm,
    )
    return results, fig


def run_wn_comparison(wn, align_shift_nm, k_linked_guess_cm,
                       k_fit_range_cm=(0.5, 6.0), x_start_cht=0.22, L_cutoff_cht=1.2,
                       L_cutoff_fit=1.5, lam0_guess_um=None, xr_range_rs=(0.22, 1.2),
                       fft_xr=None, fft_q_guess=None, data_dir='data/graphene_3x1',
                       save_dir=None, show=True):
    """
    One-call replacement for the per-wavenumber block that used to be copy-pasted
    15 times in fitting_pipeline.ipynb (data load+align -> CHT fit -> real-space
    hankel/1-sqrtx fit -> FFT comparison). Each wavenumber keeps its own cell in
    the notebook; only the per-wn tunable parameters (k_fit_range_cm,
    k_linked_guess_cm, lam0_guess_um, fft_q_guess, ...) need to be passed in.

    If save_dir is given, saves the three figures to
    {save_dir}/cht/{wn}_cht.png, {save_dir}/realspace/{wn}_realspace.png,
    {save_dir}/fft/{wn}_fft.png (directories created as needed). Set show=False
    to close figures immediately (e.g. for headless batch runs).

    Returns a dict with the CHT results plus the real-space hankel/1-sqrtx
    lambda_p/RMSE/AIC, suitable for assembling a comparison table across wn.
    """
    loaded = load_aligned_wn_signal(wn, align_shift_nm, data_dir=data_dir, L_cutoff_fit=L_cutoff_fit)
    df_target, x_f, sig_f = loaded['df_target'], loaded['x_f'], loaded['sig_f']

    if lam0_guess_um is None:
        lam0_guess_um = (2 * np.pi) / (k_linked_guess_cm * 10.0)
    if fft_q_guess is None:
        fft_q_guess = [2 * k_linked_guess_cm]
    if fft_xr is None:
        fft_xr = xr_range_rs

    cht_path = f'{save_dir}/cht/{wn}_cht.png' if save_dir else None
    if save_dir:
        os.makedirs(f'{save_dir}/cht', exist_ok=True)
        os.makedirs(f'{save_dir}/realspace', exist_ok=True)
        os.makedirs(f'{save_dir}/fft', exist_ok=True)

    cht_results, fig_cht = fit_and_plot_cht(
        x_f, sig_f, wn, x_start_cht=x_start_cht, L_cutoff_cht=L_cutoff_cht,
        k_fit_range_cm=k_fit_range_cm, k_linked_guess_cm=k_linked_guess_cm,
        save_path=cht_path)

    amplp = pd.DataFrame({'distance_um': x_f, f'{wn}_O3A': sig_f})
    outs, fig_rs, _ = compare_cavity_models(
        amplp, f'{wn}_O3A', xr=xr_range_rs, yc_um=1.9, fit_yc=False, edges='single',
        prefactors=('hankel', '1/sqrtx'), win=3, prom=0.01, lam0_guess=lam0_guess_um,
        ylim=(sig_f.min() * 2, sig_f.max() * 1.5), figsize=(8, 6))
    if save_dir:
        fig_rs.savefig(f'{save_dir}/realspace/{wn}_realspace.png', dpi=200, bbox_inches='tight')

    amp_raw, phase_raw = df_target['O3A'], df_target['O3P']
    window_len = min(41, len(amp_raw) if len(amp_raw) % 2 != 0 else len(amp_raw) - 1)
    amp_osc = amp_raw - savgol_filter(amp_raw, window_length=window_len, polyorder=2)
    fft_out = plot_channel_fft(df_target['distance_um'], amp_osc, phase_raw, label='O3', wn=wn,
                                xr=fft_xr, q_range=(0, 10), window='hann', pad_factor=3.0,
                                q_guess=fft_q_guess)
    # plot_channel_fft creates its own figure internally (via plt.subplots) rather
    # than drawing onto a pre-made one, so grab whatever it just made as "current".
    fig_fft = plt.gcf()
    if save_dir:
        fig_fft.savefig(f'{save_dir}/fft/{wn}_fft.png', dpi=200, bbox_inches='tight')

    if not show:
        plt.close(fig_cht); plt.close(fig_rs); plt.close(fig_fft)

    results = dict(wn=wn, **{k: v for k, v in cht_results.items() if k != 'wn'})
    for pf in ('hankel', '1/sqrtx'):
        p, d, met = outs[pf]['params'], outs[pf]['derived'], outs[pf]['metrics']
        damp_key = 'q_imag_um^-1' if pf == 'hankel' else 'alpha_env_um^-1'
        results[f'{pf}_lambda_p_nm'] = p['lambda_p_um'] * 1000
        results[f'{pf}_q_p_1e5cm-1'] = d['q_cm^-1'] / 1e5
        results[f'{pf}_damping'] = d['q_rad_per_um'] / p[damp_key] if p[damp_key] > 1e-9 else np.inf
        results[f'{pf}_rmse'] = met['rmse']
        results[f'{pf}_aic'] = met['aic']

    # FFT: take the dominant complex-FFT peak nearest the first q_guess as the q_p estimate
    # (matches the "q, not 2q" convention used for fft_q_guess[0] / CHT's edge-launched q).
    fft_peaks = [p for p in fft_out['peaks_complex']['peaks'] if p is not None]
    if fft_peaks:
        q_fft_1e5cm1 = fft_peaks[0] / 10.0
        results['fft_q_p_1e5cm-1'] = q_fft_1e5cm1
        results['fft_lambda_p_nm'] = (2 * np.pi / q_fft_1e5cm1) * 100
    else:
        results['fft_q_p_1e5cm-1'] = None
        results['fft_lambda_p_nm'] = None
    results['fft_damping'] = None  # not estimated from the FFT peak

    return results, (fig_cht, fig_rs, fig_fft)
