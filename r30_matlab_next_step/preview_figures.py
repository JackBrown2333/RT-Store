"""Preview the same saved data with Matplotlib; these are not native MATLAB FIGs."""
from pathlib import Path
import argparse, json
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties,fontManager


def main():
    ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('out');ap.add_argument('--font',required=True);ap.add_argument('--numeric-dir');a=ap.parse_args()
    fontManager.addfont(a.font);cn=FontProperties(fname=a.font,size=10.5)
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    m=loadmat(a.source,simplify_cells=True)['M'];t=m['time_s']-m['window_s'][0];n=len(t);dt=m['forceDt_s'];ids=m['toothIds']
    f=m['force_N'];normal=m['normalVector_N'];tc=m['moment_Nmm'][:,2]/1000
    norm=lambda x:np.linalg.norm(x,axis=-1)
    g=int(np.flatnonzero(ids==17)[0]);p=[g,int(np.flatnonzero(ids==18)[0])]
    f1=m['toothForce_N'][:,p[0],:];f2=m['toothForce_N'][:,p[1],:];pair=f1+f2
    upper=norm(f1)+norm(f2);eta=np.full(n,np.nan);mask=upper>1e-12;eta[mask]=1-norm(pair[mask])/upper[mask]
    raw=np.column_stack([f,tc,m['toothForce_N'][:,g,:],pair]);w=np.hanning(n)
    amp=2*abs(np.fft.rfft((raw-raw.mean(axis=0))*w[:,None],axis=0))/w.sum();amp[0]/=2
    if n%2==0:amp[-1]/=2
    hz=np.fft.rfftfreq(n,dt);fm=abs(m['parameters']['installation']['omegaG'])/(2*np.pi)*m['parameters']['p']['gear']['teeth']
    plt.rcParams.update({'font.family':'serif','font.serif':['Liberation Serif','DejaVu Serif'],'mathtext.fontset':'stix','font.size':10.5,
                         'legend.fontsize':9.5,'axes.linewidth':.8,'xtick.direction':'in','ytick.direction':'in',
                         'lines.linewidth':.65,'savefig.dpi':260,'axes.unicode_minus':False,
                         'axes.prop_cycle':plt.cycler(color=[(0,.145,1),(1,.11,0),(0,0,0),(.05,.55,.15)])})
    paths=[]
    captions={
      'F01_global_force':['全局力分量','全力与纯法向合力'],
      'F02_global_FFT':['力分量幅值谱','C轴转矩幅值谱'],
      'F03_tooth_17_force':['第17齿力分量','全力与纯法向合力','两齿面受力','两齿面转矩叠加'],
      'F04_pair_17_18_superposition':['单齿及双齿合力','标量上界与矢量合力','双齿合力分量','方向相消指标'],
      'F05_per_tooth_normal_contributions':['各齿纯法向合力','多齿法向载荷叠加'],
      'F06_removing_counts_and_torque':['同时参与去除的齿与齿面','C轴载荷转矩'],
      'F07_torque_time_frequency':['时域转矩','频域转矩']}
    yl={'Force (N)':'力(N)','Resultant (N)':'合力(N)','Amplitude (N)':'幅值(N)',
        'Amplitude (N m)':'幅值(N·m)','C-axis load torque (N m)':'转矩(N·m)',
        'Flank resultant (N)':'力(N)','Cancellation index':'相消指标',
        'Tooth index':'齿号','Active-removal count':'数量'}
    legend_labels={
      'F01_global_force':[['$F_x$','$F_y$','$F_z$'],['全力合力','纯法向合力']],
      'F03_tooth_17_force':[['$F_x$','$F_y$','$F_z$'],['全力合力','纯法向合力'],['齿面1','齿面2'],['齿面1','齿面2','合计']],
      'F04_pair_17_18_superposition':[['第17齿','第18齿','矢量合力'],['模长之和','矢量合力'],['$F_x$','$F_y$','$F_z$'],[]],
      'F05_per_tooth_normal_contributions':[[],['各齿模长之和','多齿矢量合力']],
      'F06_removing_counts_and_torque':[['去除齿数','去除齿面数'],[]]}
    def finish(fig,axes,name,spectrum=False):
        axes=np.asarray(axes).ravel();rows=2 if len(axes)==4 or name.startswith('F05') else 1;cols=1 if name.startswith('F05') else 2
        fig.set_layout_engine(None);fig.set_size_inches(16/2.54,(12.8 if rows==2 else 6.8)/2.54)
        colors=[(0,.145,1),(1,.11,0),(0,0,0),(.05,.55,.15)]
        for j,ax in enumerate(axes):
            row,col=divmod(j,cols);x,width=(.13,.73) if cols==1 else (.09+col*.51,.35)
            y,hh=(.30,.64) if rows==1 else (.655-row*.45,.30);ax.set_position([x,y,width,hh])
            spec=spectrum or (name.startswith('F07') and j==1)
            ax.set_xlim((0,500) if spec else (0,n*dt));ax.tick_params(direction='in',top=True,right=True,width=.8,length=3)
            ax.set_xlabel('频率(Hz)' if spec else '时间(s)',fontproperties=cn,labelpad=1)
            ax.set_ylabel(yl.get(ax.get_ylabel(),ax.get_ylabel()),fontproperties=cn,labelpad=1)
            ax.set_title('');ax.grid(False)
            for k,line in enumerate(ax.get_lines()):
                xx=np.asarray(line.get_xdata())
                if len(xx)==2 and xx[0]==xx[1]:line.set(color='.55',linestyle=':',linewidth=.6)
                else:line.set(color=colors[k%len(colors)],linewidth=.65)
            labels=legend_labels.get(name,[[] for _ in axes])[j]
            oldlegend=ax.get_legend()
            if not labels and oldlegend:labels=[z.get_text() for z in oldlegend.texts]
            if oldlegend:oldlegend.remove()
            if labels:
                components=labels==['$F_x$','$F_y$','$F_z$']
                ax.legend(ax.get_lines()[:len(labels)],labels,loc='upper right',frameon=False,prop=FontProperties(fname=a.font,size=9),handlelength=1.2 if components else 2.1,handletextpad=.35,borderpad=.1,labelspacing=.2,ncol=3 if components else 1,columnspacing=.7)
                low,high=ax.get_ylim();pad=.20 if components else .55 if len(labels)>2 else .35;ax.set_ylim(low,high+pad*(high-low))
            cy=.075 if rows==1 else y-.095
            fig.text(x+width/2,cy,'('+chr(97+j)+')'+captions[name][j],ha='center',va='center',fontproperties=cn)
        if name.startswith('F05'):
            cax=fig.axes[-1];cax.set_position([.89,.655,.02,.30]);cax.set_ylabel('力(N)',fontproperties=cn,labelpad=1)
        path=out/(name+'.png');fig.savefig(path,facecolor='white');plt.close(fig);paths.append(path)
    fig,ax=plt.subplots(2,1,figsize=(6.4,4.8),layout='constrained')
    ax[0].plot(t,f);ax[0].set(ylabel='Force (N)',title='(a) Global force components');ax[0].legend(['$F_x$','$F_y$','$F_z$'],frameon=False)
    ax[1].plot(t,np.column_stack([norm(f),norm(normal)]));ax[1].set(xlabel='Time (s)',ylabel='Resultant (N)',title='(b) Full and pure-normal vector resultants');ax[1].legend([r'$\|\sum F_i\|$',r'$\|\sum N_i\|$'],frameon=False)
    finish(fig,ax,'F01_global_force')
    fig,ax=plt.subplots(2,1,figsize=(6.4,4.8),layout='constrained')
    ax[0].plot(hz[1:],amp[1:,:3]);ax[0].set(ylabel='Amplitude (N)',title='(a) Global force: single-sided amplitude');ax[0].legend(['$F_x$','$F_y$','$F_z$'],frameon=False)
    ax[1].plot(hz[1:],amp[1:,3]);ax[1].set(xlabel='Frequency (Hz)',ylabel='Amplitude (N m)',title=f'(b) C-axis load torque; Hann, df = {hz[1]:.3f} Hz')
    for axis in ax:axis.axvline(fm,color='.5',ls='--',lw=.7)
    finish(fig,ax,'F02_global_FFT',True)
    faces=np.flatnonzero(m['faces'][:,0]==17);faceF=norm(m['faceForce_N'][:,faces,:]);faceTc=m['faceMoment_Nmm'][:,faces,2]/1000
    fig,ax=plt.subplots(2,2,figsize=(7,5.7),layout='constrained')
    ax[0,0].plot(t,m['toothForce_N'][:,g,:]);ax[0,0].set(ylabel='Force (N)',title='(a) Tooth 17: force components');ax[0,0].legend(['$F_x$','$F_y$','$F_z$'],frameon=False)
    ax[0,1].plot(t,np.column_stack([norm(m['toothForce_N'][:,g,:]),norm(m['toothNormalVector_N'][:,g,:])]));ax[0,1].set(ylabel='Resultant (N)',title='(b) Full and pure-normal resultants');ax[0,1].legend([r'$\|F\|$',r'$\|N\|$'],frameon=False)
    ax[1,0].plot(t,faceF);ax[1,0].set(xlabel='Time (s)',ylabel='Flank resultant (N)',title='(c) Individual flanks');ax[1,0].legend(['Flank 1','Flank 2'],frameon=False)
    ax[1,1].plot(t,np.column_stack([faceTc,faceTc.sum(axis=1)]));ax[1,1].set(xlabel='Time (s)',ylabel='C-axis load torque (N m)',title='(d) Flank torque contributions');ax[1,1].legend(['Flank 1','Flank 2','Sum'],frameon=False)
    finish(fig,ax,'F03_tooth_17_force')
    fig,ax=plt.subplots(2,2,figsize=(7,5.7),layout='constrained')
    ax[0,0].plot(t,np.column_stack([norm(f1),norm(f2),norm(pair)]));ax[0,0].set(ylabel='Force (N)',title='(a) Individual and combined resultants');ax[0,0].legend(['Tooth 17','Tooth 18',r'$\|F_1+F_2\|$'],frameon=False)
    ax[0,1].plot(t,np.column_stack([upper,norm(pair)]));ax[0,1].set(ylabel='Force (N)',title='(b) Scalar upper bound and vector sum');ax[0,1].legend([r'$\|F_1\|+\|F_2\|$',r'$\|F_1+F_2\|$'],frameon=False)
    ax[1,0].plot(t,pair);ax[1,0].set(xlabel='Time (s)',ylabel='Force (N)',title='(c) Combined force components');ax[1,0].legend(['$F_x$','$F_y$','$F_z$'],frameon=False)
    ax[1,1].plot(t,eta,color='#cc79a7');ax[1,1].set(xlabel='Time (s)',ylabel='Cancellation index',ylim=(0,1),title='(d) Directional cancellation index')
    finish(fig,ax,'F04_pair_17_18_superposition')
    perN=norm(m['toothNormalVector_N'])
    fig,ax=plt.subplots(2,1,figsize=(6.4,4.8),layout='constrained')
    im=ax[0].imshow(perN.T,origin='lower',aspect='auto',extent=(t[0],t[-1],ids[0]-.5,ids[-1]+.5),cmap='viridis');ax[0].set(ylabel='Tooth index',yticks=ids,title='(a) Saved per-tooth pure-normal resultants');fig.colorbar(im,ax=ax[0],label='Force (N)')
    ax[1].plot(t,np.column_stack([perN.sum(axis=1),norm(normal)]));ax[1].set(xlabel='Time (s)',ylabel='Force (N)',title='(b) Multitooth vector cancellation');ax[1].legend([r'$\sum\|N_i\|$',r'$\|\sum N_i\|$'],frameon=False)
    finish(fig,ax,'F05_per_tooth_normal_contributions')
    fig,ax=plt.subplots(2,1,figsize=(6.4,4.8),layout='constrained')
    ax[0].step(t,m['nRemovingTeeth'],where='post');ax[0].step(t,m['nRemovingFaces'],where='post');ax[0].set(ylabel='Active-removal count',ylim=(0,7),title='(a) Simultaneous material removal');ax[0].legend(['Teeth','Flanks'],frameon=False)
    ax[1].plot(t,tc);ax[1].set(xlabel='Time (s)',ylabel='C-axis load torque (N m)',title='(b) Signed load torque; fixed rotation direction')
    finish(fig,ax,'F06_removing_counts_and_torque')
    fig,ax=plt.subplots(1,2);ax[0].plot(t,tc);ax[0].set_ylabel('C-axis load torque (N m)');ax[1].plot(hz[1:],amp[1:,3]);ax[1].set_ylabel('Amplitude (N m)');ax[1].axvline(fm,color='.55',ls=':')
    finish(fig,ax,'F07_torque_time_frequency')
    fig=plt.figure(figsize=(13,10.5));places=[(.035,.735,.93,.245),(.035,.475,.93,.245),(.035,.015,.45,.445),(.515,.015,.45,.445)]
    for pos,path in zip(places,[paths[0],paths[6],paths[2],paths[3]]):
        axis=fig.add_axes(pos);axis.imshow(plt.imread(path));axis.axis('off')
    fig.savefig(out/'r30_academic_preview.png',dpi=190,facecolor='white');plt.close(fig)
    # Cross-check the independent previews against the actual MATLAB-compatible export.
    matlab=Path(a.numeric_dir)/'fft_single_sided_amplitude.csv' if a.numeric_dir else Path('no_numerical_reference')
    if matlab.exists():
        saved=np.loadtxt(matlab,delimiter=',',skiprows=1)
        error=float(np.max(abs(saved[:,1:]-amp)))
        assert error<1e-10
        (out/'preview_validation.json').write_text(json.dumps({'renderer':'Matplotlib, not MATLAB',
          'same_source_M':True,'FFT_max_difference_from_Octave_export':error,
          'pair_peak_N':float(norm(pair).max()),'scalar_sum_peak_N':float(upper.max()),
          'n_samples':n,'sampleRate_Hz':1/dt,'frequency_bin_Hz':hz[1],
          'style_reference_pdf_pages':[40,54],'plot_time_origin_s':float(m['window_s'][0]),
          'Chinese_preview_font':'Noto Sans CJK SC fallback; native MATLAB requests SimSun'},indent=2))

if __name__=='__main__':main()
