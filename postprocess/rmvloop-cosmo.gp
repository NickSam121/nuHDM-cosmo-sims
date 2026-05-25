reset
set encoding utf8
gifdelay=10
lws=2.
flgjpg=0
rgbborder='#000000'
rgblabel='#00ff00'
outtrunc='rmv'
vardens=1
varscale=1
vartextcolor=1
##### TEXTCUBEHELIX SETTINGS: ####
## Try to set phi_start to estimated inverse colorscale end value and rotations reversed
## to maximise contrast.
## Colorhelix to inverse cubehelix:
## mod(deg_end-90deg, 360) / 120 = mod(deg_start+rot*360, 360) / 120
## OR use flginvert=1. Manual adjustment may be useful.
phi_start=2.    # start color (1.0=red, 2.0=green etc.; modulo 3)
rotations=0.    # number of colour space rotations along grayscale
sathue=2.0      # saturation of hue (0=gray...1=colour)
gamma_gry=1.0   # gamma correction

#### Gray range
x0=0.; x1=1.    # range of grayscale (standard = 0...1)
#### Non-standard exponents (not used by Green (2011)
gamma_rot=1.    # "gamma" correction for rotation rate
power_amp=1.    # optional exponent for colour amplification (<0 for lin. trunc.)

#### Expert-options
sw_sc=-1.       #s-curve parameter.
power_sc=8      #
kr=1.; kg=1.; kb=1.; or=0.; og=0.; ob=0.
flgdistort=0    # experimental distortion mode
flginvert=0     # inverse colours
#### RGB weights ####
#wr=0.30; wg=0.59; wb=0.11  # NTSC
#wr=0.299; wg=0.587; wb=0.114 # PAL
wr=0.2126;wg=0.7152;wb=0.0722# HDTV
#wr=1.;wg=1.;wb=1.
pwr_aexp=1.
load 'cubehelix_x.gp'
########

set term jpeg enhanced size 1500,1200 font arial 24 background '#000000'; rgbborder='#ffffff'; flgjpg=1
set border lw lws lc rgb rgbborder
set tics in front
#set tics out
set size ratio -1
set style line 1 pt 1 lt 1 lw 1.0*lws lc rgb '#ff0000'
set style line 2 pt 2 lt 1 lw 1.0*lws lc rgb '#00ff00'
set style line 3 pt 3 lt 1 lw 1.0*lws lc rgb '#0000ff'
set style line 4 pt 1 lt 1 lw 1.0*lws lc rgb '#00ffff'
set style line 5 pt 2 lt 1 lw 1.0*lws lc rgb '#ff00ff'
set style line 6 pt 3 lt 1 lw 1.0*lws lc rgb '#ffff00'
set style line 7 pt 3 lt 1 lw 1.0*lws lc rgb '#a0a0a0'

set style line 10 pt 1 lt 1 lw 1.0*lws lc rgb '#ffff88'

lmg=0.11
rmg=0.86
bmg=0.05
tmg=0.98

set lmargin screen lmg
set rmargin screen rmg
set bmargin screen bmg
set tmargin screen tmg

unset logscale xy
#set term x11 2 title "Cluster 3D" size 720,720
kzoom=1.6
set xyplane 0
set view equal xyz
set view 0,0,kzoom
fbeta(r,t)=1.-(t/r)**2

## Constants
kB=1.3806200e-16
mH=1.6600000e-24
## Units
u2ast=1.#cgs
#u2ast=4.7884#cgs to 10^9 Msun / kpc
u2kms=65.5926606#Nbody to km/s
xyscale=1.

## 3D positions
#w=1.e3; off=0.e3
#off=w
off=0.e0
offx=0.; offy=0.; offz=0.
w=1.; dxy=0.1; nmtics=5
hw=0.5*w

#set xrange [-hw+off:hw+off]; set yrange [-hw+off:hw+off]; set zrange [-hw+off:hw+off]
set xrange [-hw+offx:hw+offx]; set yrange [-hw+offy:hw+offy]; set zrange [-hw+offz:hw+offz]
#set xrange [0:2*w]; set yrange [0:2*w]; set zrange [0:2*w]
#set xrange [500:1500]; set yrange [500:1500]; set zrange [500:1500]
#set xlabel 'x kpc' tc rgb rgbborder; set ylabel 'y kpc' offset 1,0 tc rgb rgbborder; set zlabel 'z kpc' tc rgb rgbborder
#set xlabel 'x / Mpc h^{-1}' tc rgb rgbborder; set ylabel 'y / Mpc h^{-1}' offset 1,0 tc rgb rgbborder
#set ylabel 'x,y / Mpc' offset 1,0 tc rgb rgbborder
set label 1 'x,y / Mpc' at screen 0.03,0.47 rotate by 90 tc rgb rgbborder
#unset cblabel
#set cblabel 'log {/Symbol S}' offset 1,0 tc rgb '#7f7f7f'
set cbtics 1
set cbtics 20 #velocity
set mcbtics 5
set xtics dxy; set ytics dxy
set mxtics nmtics; set mytics nmtics

set palette file './pycolorhelix_nebula0001.rgba'

cbmin=0.; cbmax=0.1
#cmin=-300.; cmax=300.
cmin=-100.; cmax=100.
ugamma=0.25
#lgsmin=-7.; lgsmax=-6.5; set cbtics 0.1; set mcbtics 5#fullbox
#lgsmin=-0.1; lgsmax=2.0; set cbtics 0.5; set mcbtics 5#cosmo without vardens
#lgsmin=-10; lgsmax=-1.0; set cbtics 0.5; set mcbtics 5#cosmo with vardens
#swd=0.9; clg=-6.; dlg=4.; lgsmin=clg+(swd-1.)*dlg; lgsmax=clg+swd*dlg; set cbtics 0.1*dlg; set mcbtics 5
#lgsmin=1; lgsmax=7; set cbtics 0.5; set mcbtics 5#cosmo temp
#clg= 4.; dlg=0.1; lgsmin=clg-0.5*dlg; lgsmax=clg+0.5*dlg; set cbtics 0.1*dlg; set mcbtics 5 #temperature
#lgsmin=-2.; lgsmax=2. #pressure
#lgsmin=0.; lgsmax=2. #velocity
lgsmin=-17.; lgsmax=-8. #background particles
smin=10.**lgsmin; smax=10.**lgsmax
flog(x)=x<smin ? lgsmin : log10(x)
i2myr=100; i2off=0
#i2myr=100; i2off=-8600

####
rmvpath='./movie1/'
####

#### Don't change
#plotmode=1#power-law density
plotmode=2#log density
#plotmode=0#velocity

#### Select data type
#chvar='dens';kscal=u2ast; set cblabel 'log {/Symbol S}' offset 1,0 tc rgb rgbborder
#chvar='temp';kscal=1.; set cblabel 'log T' offset 1,0 tc rgb '#7f7f7f'
#chvar='pres';kscal=1.; set cblabel 'log P' offset 1,0 tc rgb '#7f7f7f'
chvar='var6';kscal=u2ast
#chvar='stars';kscal=u2ast
#chvar='vx'; kscal=65.5926606
####
if (plotmode==0) cbmin=cmin; cbmax=cmax
if (plotmode==1) cbmin=0.; cbmax=cbmax**ugamma
if (plotmode==2) cbmin=lgsmin; cbmax=lgsmax
set cbrange [cbmin:cbmax]
#if (plotmode/10==3) set cbrange [0:255]
nfile=1100
istep=1
ifile0=0
iout=-1

set key left

do for [ifile=ifile0:nfile:istep] {
    iout=iout+1
    nchar=sprintf("%05d",ifile)
    nchout=sprintf("%05d",iout)
    infofile=rmvpath.'info_'.nchar.'.txt'
    load infofile
    unit_v=unit_l/unit_t
    unit_T2=mH/kB * unit_v**2
    if (ifile==ifile0) {aexp0=aexp}
    tmyr=time*unit_t/31.5576e13
    zred=1./aexp-1.
    chaexp=sprintf("a = %.3f",aexp)
    chzred=sprintf("z = %.3f",zred)
    #chmyr=sprintf("%5d Myr",tmyr)
    chlabel=chaexp.'; '.chzred
    if (chvar eq 'dens') {unit_plot=1.e24*unit_d}
    #if (chvar eq 'temp') {unit_plot=1.e-9*unit_T2}
    if (chvar eq 'temp') {unit_plot=kscal}
    if (chvar eq 'stars') {unit_plot=kscal}
    if (vardens>0) {kscal=unit_plot}
    #if (vardens>0 && ifile==ifile0) {print sprintf('unit_plot = %.6e',unit_plot)}
    if (vardens>1) {set autoscale cb; set cbtics autofreq}
    if (varscale>0) {xyscale=unit_l/3.085677581282e24; hw=0.5*xyscale
	set xrange[-hw:hw]; set yrange[-hw:hw]
	set xtics autofreq; set ytics autofreq; unset mxtics; unset mytics
    }
    rmvfile=rmvpath.chvar.'_'.nchar.'.map'
    #rmvfile2=rmvpath.chvar2.'_'.nchar.'.map'
    jpgoutput=outtrunc.nchout.'.jpg'
    #ismissing=system("ismissing.sh ".datafile)
    file_exists(file) = system("[ -f '".file."' ] && echo '1' || echo '0'") + 0
    ismissing=1-file_exists(rmvfile)
    jpgexists=file_exists(jpgoutput)
    if (jpgexists) {print 'i,jpgoutput: ',ifile,' ',jpgoutput,' exists'; continue}
    print 'i,rmv,ismissing =',iout,' ',rmvfile,' ',ismissing
    if (vartextcolor>0) {
	if (flginvert>0) {
	    xgray=1.-(aexp-aexp0)**pwr_aexp
	} else {
	    xgray=(aexp-aexp0)**pwr_aexp
	}
	rgblabel=f2rgbhex(xgray)
	#print('xgray    = ',xgray)
	#print('rgblabel = ',rgblabel)
	#print sprintf("xgray = %.3f, rgblabel = %s",xgray,rgblabel)
    }
    set label 2 chlabel at graph 0.05,0.96,0.96 front tc rgb rgblabel
    if (ismissing) {unset output; q}
    if (flgjpg>0) {set output jpgoutput}
    if (plotmode==0||plotmode==3) {plot rmvfile u ($1*xyscale):($2*xyscale):($3*kscal) w image}
    if (plotmode==1) {plot rmvfile u ($1*xyscale):($2*xyscale):(($3*kscal)**ugamma) w image}
    if (plotmode==2) {plot rmvfile u ($1*xyscale):($2*xyscale):(flog($3*kscal)) w image}
}
