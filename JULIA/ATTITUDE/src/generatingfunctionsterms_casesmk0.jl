######
# terms of the generating functions in the special case m=0 used in generatingfunctionsrpdrag_p1.jl; generatingfunctiongravity.jl generatingfunctionmagnetic.jl
#####

function getgeneratinfunctiontermsW2_allpert_psil_psig_smk0(psil,psig,npsil,npsig)
    VG = zeros(8);
    VGG = zeros(8);
    psigl = psil+psig;

    t1 = sin(psigl)
    t3 = 2 * npsil + 2 * npsig
    t4 = 1 / t3
    VG[1] = t4 * t1
    t5 = 2 * psigl
    t6 = sin(t5)
    VG[2] = t4 * t6 / 2
    t8 = 3 * psigl
    t9 = sin(t8)
    VG[3] = t4 * t9 / 3
    t11 = 4 * psigl
    t12 = sin(t11)
    VG[4] = t4 * t12 / 4
    t14 = cos(psigl)
    VG[5] = -t4 * t14
    t16 = cos(t5)
    VG[6] = -t4 * t16 / 2
    t19 = cos(t8)
    VG[7] = -t4 * t19 / 3
    t22 = cos(t11)
    VG[8] = -t4 * t22 / 4
    t25 = t3 ^ 2
    t26 = 1 / t25
    VGG[1] = -t14 * t26
    VGG[2] = -t16 * t26 / 4
    VGG[3] = -t19 * t26 / 9
    VGG[4] = -t22 * t26 / 16
    VGG[5] = -t1 * t26
    VGG[6] = -t6 * t26 / 4
    VGG[7] = -t9 * t26 / 9
    VGG[8] = -t12 * t26 / 16
    VG  = 2.0*VG;
    VGG = 4.0*VGG;
    return VG,VGG;
    
end