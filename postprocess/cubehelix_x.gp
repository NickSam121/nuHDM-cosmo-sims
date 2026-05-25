########################################################################
##                 **** CUBE HELIX COLOUR SCHEME ****                 ##
##                      for GNUPLOT 5.4 or later                      ##
##           Based on a method described by D.A. Green (2011)         ##
##           http://arxiv.org/abs/1108.5083                           ##
##                                                                    ##
## Modified by Ingo Thies in 2022 for automatic text color adjustment ##
##                                                                    ##
## Note: cubehelix.gp is not yet compatible with colormanager.gp,     ##
## but the output cubehelix.rgb can be used with it.                  ##
########################################################################
##### EXAMPLE SETTINGS: ####
#phi_start=0.5   # start color (1.0=red, 2.0=green etc.; modulo 3)
#rotations=-1.5  # number of colour space rotations along grayscale
#sathue=1.2      # saturation of hue (0=gray...1=colour)
#gamma_gry=1.0   # gamma correction
##### Non-standard exponents (not used by Green (2011)
#gamma_rot=1.    # "gamma" correction for rotation rate
#power_amp=1.    # optional exponent for colour amplification (<0 for lin. trunc.).
##### Expert options
#sw_sc=0.        # s-curve parameter.
#power_sc=1.     #
#kr=1.; kg=1.; kb=1.; or=0.; og=0.; ob=0.
#flgdistort=0    # experimental distortion mode
flginvert=0      # inverse colours
##### RGB weights ####
#wr=0.30; wg=0.59; wb=0.11  # NTSC
#wr=0.299; wg=0.587; wb=0.114 # PAL
#wr=0.2126;wg=0.7152;wb=0.0722# HDTV
#wr=1.;wg=1.;wb=1
########
#### Some extra definitions
fcrop01(x)= x>1. ? 1 : x<0. ? 0. : x
########
wsum=wr+wg+wb
wr=wr/wsum;wg=wg/wsum;wb=wb/wsum
#### u,v vectors used by Green (2011)
#u1=-0.074306;u2=-0.146135;u3=0.986470
#v1=0.891385; v2=-0.453247;v3=0.
######## Calculate the colour unit vectors here
sqrho=wr*wr+wg*wg
rhorg=sqrt(sqrho)
eps0=wb*wr/sqrho
unorm=sqrt(1.+eps0**2+(wg/wr*eps0)**2)
epsu=eps0/unorm
u1=-epsu;u2=-epsu*wg/wr;u3=1./unorm
v1=wg/rhorg;v2=-wr/rhorg;v3=0
#### Expert mode: stretch & shift base vectors
ku1=kr;ku2=kg;ku3=kb; kv1=ku1;kv2=ku2;kv3=ku3
u1off=or;u2off=og;u3off=ob; v1off=u1off;v2off=u2off;v3off=u3off
if (flgdistort>=1) {
    print 'Expert mode: stretch & shift enabled'
    u1=u1*ku1+u1off;u2=u2*ku2+u2off;u3=u3*ku3+u3off
    v1=v1*kv1+v1off;v2=v2*kv2+v2off;v3=v3*kv3+v3off
}
#print 'u1,u2,u3 = ',sprintf("%+.6f, %+.6f, %+.6f",u1,u2,u3)
#print 'v1,v2,v3 = ',sprintf("%+.6f, %+.6f, %+.6f",v1,v2,v3)
########
#### Function definitions (note that the meaning of 'x' varies
#    to reduce number of function calls)
## s-curve mode
fx05(x)=2.*(x-0.5)
if (sw_sc<0 && power_sc>=0) {
    fxs(x)=(1.-abs(fx05(x)))**power_sc
    fxw(x)=1.+sgn(x-0.5)*(1.-fxs(x))
    fsc(x)=0.5*abs(sw_sc)*fxw(x)+(1.-abs(sw_sc))*x
} else {
    fxs(x)=abs(fx05(x))**abs(power_sc) * sgn(fx05(x))
    fxw(x)=sw_sc*fxs(x)+(1.-sw_sc)*fx05(x)
    fsc(x)=0.5*(fxw(x)+1.)
}
## Gray, rot & amp
fxgry(x)=fsc(x)**gamma_gry
fxrot(x)=x**gamma_rot
#fxamp(x)=x**power_amp
fphi(x)=2.*pi*(phi_start/3.+1.+rotations*x)
dev(x)=x*(1.-x)
pdev(x)=0.25*(4*dev(x))**power_amp

## the following functions is an experimental linear scaling for more saturated colors
rdev(x)=x<0.5 ? x : 1.-x
prdev(x)=0.4*(2.5*rdev(x))**(-power_amp)
##
famp(x)=sathue*(power_amp==1 ? dev(x) : power_amp>0 ? pdev(x) : prdev(x))
fuv1(x)=u1*cos(x)+v1*sin(x)
fuv2(x)=u2*cos(x)+v2*sin(x)
fuv3(x)=u3*cos(x)+v3*sin(x)
f01(x,y)=x+famp(x)*fuv1(y)
f02(x,y)=x+famp(x)*fuv2(y)
f03(x,y)=x+famp(x)*fuv3(y)
f0red(x)=f01(fxgry(x),fphi(fxrot(x)))
f0grn(x)=f02(fxgry(x),fphi(fxrot(x)))
f0blu(x)=f03(fxgry(x),fphi(fxrot(x)))
if (flginvert>0) {
    fred(x)=1.-fcrop01(f0red(x))
    fgrn(x)=1.-fcrop01(f0grn(x))
    fblu(x)=1.-fcrop01(f0blu(x))
} else {
    fred(x)=fcrop01(f0red(x))
    fgrn(x)=fcrop01(f0grn(x))
    fblu(x)=fcrop01(f0blu(x))
}
fgry(x)=wr*fred(x)+wg*fgrn(x)+wb*fblu(x)
#### Optional linecolor functions (requires gnuplot 4.5+)
frgb(x)=int(255*fred(x))*65536 + int(255*fgrn(x))*256 + int(255*fblu(x))
fggg(x)=int(255*fgry(x))*65536 + int(255*fgry(x))*256 + int(255*fgry(x))
rgbhex(r,g,b)=sprintf("#%06x",65536*r+256*g+b)        #0 to 255
#f2rgbhex(x)=rgbhex(fred(x),fgrn(x),fblu(x))
f2rgbhex(x)=sprintf("#%06x",frgb(x))
