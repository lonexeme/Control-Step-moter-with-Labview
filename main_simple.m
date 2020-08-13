% 7段S型加减速（一组数据）
clc;
clear;
close all;
% 输入参数
q0 = 0; q1 = 10;
v0 = 0; v1 = 0;
vmax = 10; amax = 10; jmax = 30;
% 计算个阶段参数Ta, Tv, Td, Tj1, Tj2, q_0, q_1, v_0, v_1, vlim, j_max, j_min
calresult = calparas(q0,q1,v0,v1,vmax,amax,jmax);
% 运行总时间
T = calresult(1) + calresult(2) + calresult(3);
i = 1;
% 计算T时间内每1ms的p、vel、acc、jerk参数
for t = 0: 0.001: T
    time(i) = 0.001*i;
    data_matrix(i,:) = caloutdata(t, calresult(1), calresult(2), calresult(3), calresult(4), calresult(5), calresult(6), calresult(7), ...
                       calresult(8), calresult(9), calresult(10), calresult(11), calresult(12), calresult(13));
    i = i + 1;
end

time=time*1000;
for m = 1: 4
    subplot(4,1,m)
    plot(time, data_matrix(:,m), 'LineWidth', 2);
    axis tight
    grid on
end




