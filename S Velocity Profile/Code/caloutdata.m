function dataresult = caloutdata(t, Ta, Tv, Td, Tj1, Tj2, q0, q1, v0, v1, vlim, jmax, jmin, sigma)
    % 计算当前时间段位置、速度、加速度、加加速度值
    % 此处显示详细说明
    % Ta加速阶段持续时间
    % Tv匀速阶段持续时间
    % Td减速阶段持续时间
    % Tj1加加速阶段持续时间
    % Tj2减加速阶段持续时间
    % q_0, q_1, v_0, v_1, a_max, a_min, j_max, j_min均为转换后的参数
    % vlim为当前阶段所能达到的最大速度
        
    % 运行总时间
    T = Ta + Tv + Td;
    % 加速段
    if (t >= 0 && t < Tj1)
        q = q0 + v0*t + jmax*power(t,3)/6;
        q_dev = v0 + jmax*(t^2/2);
        q_dev2 = jmax*t;
        q_dev3 = jmax;
    elseif (t >= Tj1 && t < Ta - Tj1)
        alima = jmax*Tj1;
        q = q0 + v0*t +(alima/6)*(3*t^2 - 3*Tj1*t + Tj1^2);        
        q_dev = v0 + alima*(t - Tj1/2);
        q_dev2 = alima;
        q_dev3 = 0;
    elseif (t >= Ta - Tj1 && t < Ta)
        q = q0 + (vlim + v0)*(Ta/2) - vlim*(Ta - t) - jmin*(power(Ta - t,3)/6);
        q_dev = vlim + jmin*(power(Ta - t, 2)/2);
        q_dev2 = -jmin*(Ta - t);
        q_dev3 = jmin;
        % 匀速段
    elseif (t >= Ta && t < Ta + Tv)
        q = q0 + (vlim + v0)*(Ta/2) + vlim*(t - Ta);
        q_dev = vlim;
        q_dev2 = 0;
        q_dev3 = 0;
        % 减速段
    elseif (t >= Ta + Tv && t < T - Td + Tj2)
        q = q1 - (vlim + v1)*(Td/2) + vlim*(t - T + Td) - jmax*(power(t - T + Td, 3)/6);
        q_dev = vlim - jmax*(power(t - T + Td, 2)/2);
        q_dev2 = -jmax*(t - T + Td);
        q_dev3 = -jmax;
    elseif (t >= T - Td + Tj2 && t < T - Tj2)
        alimd = -jmax*Tj2;
        q = q1 - (vlim + v1)*(Td/2) + vlim*(t - T + Td) + (alimd/6)*(3*power(t - T + Td, 2) - 3*Tj2*(t - T + Td) + Tj2^2);
        q_dev = vlim + alimd*(t - T + Td - Tj2/2);
        q_dev2 = alimd;
        q_dev3 = 0;
    elseif (t >= T - Tj2 && t <= T)
        q = q1 - v1*(T - t) - jmax*(power(T - t, 3)/6);
        q_dev = v1 + jmax*(power(t - T, 2)/2);
        q_dev2 = -jmax*(T - t);
        q_dev3 = jmax;
    end
    
    % 按照公式（3.33）进行输出转换
    
    q_out = sigma*q;
    q_dev_out = sigma*q_dev;
    q_dev2_out = sigma*q_dev2;
    q_dev3_out = sigma*q_dev3;
    dataresult = [q_out,q_dev_out,q_dev2_out,q_dev3_out];
    
end

