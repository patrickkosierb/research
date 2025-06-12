
init

%% no noise+ states known
run_ideal

tiledlayout(2,1);
nexttile;
plot(t,j1_ext)
title("ext1")
legend("x","y","z")
nexttile;
plot(t,tau_ext_ideal(1,:))
title("int1")
legend("z")

figure;
tiledlayout(2,1);
nexttile;
plot(t,j2_ext)
title("ext2")
legend("x","y","z")
nexttile;
plot(t,tau_ext_ideal(2,:))
title("int2")
legend("z")

figure;
tiledlayout(2,1);
nexttile;
plot(t,j3_ext)
title("ext3")
legend("x","y","z")
nexttile;
plot(t,tau_ext_ideal(3,:))
title("int3")
legend("z")

figure;
tiledlayout(2,1);
nexttile;
plot(t,j4_ext)
title("ext4")
legend("x","y","z")
nexttile;
plot(t,tau_ext_ideal(4,:))
title("int4")
legend("z")

figure;
tiledlayout(2,1);
nexttile;
plot(t,j5_ext)
title("ext5")
legend("x","y","z")
nexttile;
plot(t,tau_ext_ideal(5,:))
title("int5")
legend("z")

figure;
tiledlayout(2,1);
nexttile;
plot(t,j6_ext)
title("ext6")
legend("x","y","z")
nexttile;
plot(t,tau_ext_ideal(6,:))
title("int6")
legend("z")


figure;
plot(t,q_true)