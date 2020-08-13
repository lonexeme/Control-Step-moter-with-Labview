% 7段S型加减速（一组数据）
clc;
clear;
close all;
% 输入参数
data = [0, 10, 20, 60, 80, -10, 20, -30, -100];

% 数据长度
len = length(data);
dir = sign(diff(data));

% 设置前瞻速度，也可以自定义速度
setvel = 5;
for i=1:len
    if i==1
        q0(i)=data(i);
        q1(i)=data(i+1);
        v0(i)=0;
        if dir(i)<0
            v1(i)=-setvel;
        else
            v1(i)=setvel;
        end
    elseif i==len
        v1(i-1)=0;
    else
        q0(i)=data(i);
        q1(i)=data(i+1);
        v0(i)=v1(i-1);
        if dir(i)<0
            v1(i)=-setvel;
        else
            v1(i)=setvel;
        end
    end
end

% 设置输入参数
vmax = 100; amax = 1000; jmax = 3000;
for j =1:len-1
    % 计算个阶段参数Ta, Tv, Td, Tj1, Tj2, q_0, q_1, v_0, v_1, vlim, j_max, j_min
    calresult(j,:) = calparas(q0(j), q1(j), v0(j), v1(j), vmax, amax, jmax);
    % 运行总时间
    T(j) = calresult(j,1) + calresult(j,2) + calresult(j,3);
end

% 定义变量n 用于连接数据
n = 0;
for k = 1:len-1
    m = 1;
    % 计算T时间内每1ms的p、vel、acc、jerk参数，并存放于data_matrix矩阵中
    for t = 0: 0.001: T(k)
        time(m) = 0.001*m;
        data_matrix(m,:) = caloutdata(t, calresult(k,1), calresult(k,2), calresult(k,3), calresult(k,4), calresult(k,5), calresult(k,6), calresult(k,7), ...
                           calresult(k,8), calresult(k,9), calresult(k,10), calresult(k,11), calresult(k,12), calresult(k,13));
        m = m + 1;
    end
    n = n + m;
    if k == 1
        % 创建时间数组
        time_one(1:m-1) = time;
        % 记录时间节点，用于连接时间数组
        time_end = time(end);
        % 分别从data_matrix矩阵中读取位置、速度、加速度、加加速度数据
        data_out_pos(1:m-1) = data_matrix(:,1)';
        data_out_vel(1:m-1) = data_matrix(:,2)';
        data_out_acc(1:m-1) = data_matrix(:,3)';
        data_out_jerk(1:m-1) = data_matrix(:,4)';
        start_index = m;
        % 清除工作区变量，避免由于矩阵不匹配引起的报错
        clear time data_matrix
    else
        time_one(start_index:n-k) = time + time_end;
        time_end = time_one(end);
        data_out_pos(start_index:n-k) = data_matrix(:,1)';
        data_out_vel(start_index:n-k) = data_matrix(:,2)';
        data_out_acc(start_index:n-k) = data_matrix(:,3)';
        data_out_jerk(start_index:n-k) = data_matrix(:,4)';
        start_index = n - k + 1;
        clear time data_matrix
    end
end

% 将多维矩阵转换为一维矩阵
data_matrix_one = [data_out_pos',data_out_vel',data_out_acc',data_out_jerk'];
time_one = time_one*1000;

% plot输出图标数据
for h = 1: 4
    subplot(4,1,h)
    plot(time_one, data_matrix_one(:,h), 'LineWidth', 2);
    axis tight
    grid on
end




