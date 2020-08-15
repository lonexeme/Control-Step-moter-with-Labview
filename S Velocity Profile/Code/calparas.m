function calresult = calparas(q0,q1,v0,v1,vmax,amax,jmax)
    % 由输入参数计算所需的各项参数
    % 输入参数
    % q0初始位置参数
    % q1终点位置参数（可作为下一段的初始位置参数使用）
    % v0初始速度参数
    % v1终点速度参数
    % vmax运行阶段最大能达到的速度
    % amax运行阶段最大能达到的加速度
    % jmax运行阶段最大能达到的加加速度
    % 输出参数
    % Ta加速阶段持续时间
    % Tv匀速阶段持续时间
    % Td减速阶段持续时间
    % Tj1加加速阶段持续时间
    % Tj2减加速阶段持续时间
    % q_0, q_1, v_0, v_1, a_max, a_min, j_max, j_min均为转换后的参数
    % vlim为当前阶段所能达到的最大速度
    
    
    %% 从3.4.3节中易知使用取相反数的方式来赋值vmin、amin、jmin
    vmin = -vmax; amin = -amax; jmin = -jmax;
    %% 结合公式（3.31）-（3.32）将输入参数转化为计算所需的实际参数
    sigma = sign(q1 - q0);
    q_0 = sigma*q0;
    q_1 = sigma*q1;
    v_0 = sigma*v0;
    v_1 = sigma*v1;
    v_max = ((sigma+1)/2)*vmax + ((sigma-1)/2)*vmin;
    v_min = ((sigma+1)/2)*vmin + ((sigma-1)/2)*vmax;
    a_max = ((sigma+1)/2)*amax + ((sigma-1)/2)*amin;
    a_min = ((sigma+1)/2)*amin + ((sigma-1)/2)*amax;
    j_max = ((sigma+1)/2)*jmax + ((sigma-1)/2)*jmin;
    j_min = ((sigma+1)/2)*jmin + ((sigma-1)/2)*jmax;
    
    
    %% 结合公式（3.21）-（3.25）计算参数Ta、Tv、Td、Tj1、Tj2
    % Case1->vlim=vmax,判断是否能够达到最大加速度
    if (v_max - v_0)*j_max < a_max^2
        % 达不到 a_max
        Tj1 = sqrt((v_max - v_0) / j_max);
        Ta = 2*Tj1;
    else
        % 达到a_max
        Tj1 = a_max / j_max;
        Ta = Tj1 + (v_max - v_0) / a_max;
    end
    if (v_max - v_1)*j_max < a_max^2
        % 达不到a_min
        Tj2 = sqrt((v_max - v_1) / j_max);
        Td = 2*Tj2;
    else
        % 达到a_min
        Tj2 = a_max / j_max;
        Td = Tj2 + (v_max - v_1) / a_max;
    end
    
    Tv = (q_1 - q_0)/v_max - (Ta/2)*(1 + v_0/v_max) - (Td/2)*(1 + v_1/v_max);
    
    % 对Tv的值进行判断
    if Tv > 0
        % 达到vmax
        vlim = v_max;
        T = Ta + Tv + Td;
        calresult = [Ta, Tv, Td, Tj1, Tj2, q_0, q_1, v_0, v_1, vlim, j_max, j_min, sigma];
    else
        % Case2->vlim<vmax,不存在匀速状态，即Tv=0
        % 结合公式（3.26）-（3.27）计算参数Ta、Td
        % 此段未避免程序复杂性，可做子函数，为后续做代码复用提供便利。此处为说明逻辑故未作处理
        Tv = 0;
        Tj = a_max / j_max;
        Tj1 = Tj;
        Tj2 = Tj;
        delta = (a_max^4/j_max^2) + 2*(v_0^2 + v_1^2) + a_max*(4*(q_1 - q_0) - 2*(a_max/j_max)*(v_0 + v_1));
        Ta = ((power(a_max, 2)/j_max) - 2.0*v_0 + sqrt(delta)) / (2.0*a_max);
        Td = ((power(a_max, 2)/j_max) - 2.0*v_1 + sqrt(delta)) / (2.0*a_max);
        
        % 对Ta和Td进行讨论
        if Ta < 0 || Td < 0
            if Ta < 0 && v_0 >= v_1
                % 只有减速段,Ta=0，结合公式（3.28）计算
                Ta = 0;
                Tj1 = 0;
                Td = 2*(q_1 - q_0) / (v_0 + v_1);
                Tj2 = (j_max*(q_1 - q_0) - sqrt(j_max*(j_max*power(q_1 - q_0, 2) + power(v_1 + v_0, 2)*(v_1 - v_0)))) / (j_max*(v_1 + v_0));
                vlim = v_0;
            elseif Td < 0 && v_1 >= v_0
                % 只有加速段，Td=0，结合公式（3.29）计算
                Td = 0;
                Tj2 = 0;
                Ta = 2*(q_1 - q_0) / (v_0 + v_1);
                Tj1 = (j_max*(q_1 - q_0) - sqrt(j_max*(j_max*power(q_1 - q_0, 2)) - power(v_1 + v_0, 2)*(v_1 - v_0))) / (j_max*(v_1 + v_0));
                a_lima = j_max * Tj1;
                vlim = v_0 + (Ta - Tj1) * a_lima;
            end
            calresult = [Ta, Tv, Td, Tj1, Tj2, q_0, q_1, v_0, v_1, vlim, j_max, j_min, sigma];
        elseif  Ta > 2*Tj && Td > 2*Tj
            % 加速和减速阶段均能达到最大加速度
            a_lima = j_max * Tj1;
            vlim = v_0 + a_lima*(Ta - Tj);
            calresult = [Ta, Tv, Td, Tj1, Tj2, q_0, q_1, v_0, v_1, vlim, j_max, j_min, sigma];
        else
            % 加速和减速阶段至少有一段不能达到最大加速度，不常见。
            % 精确值很难计算，采用近似值来逼近，虽然不是最优解，但是满足使用要求
            gamma = 0.99;
            % 0<gamma<1，越接近1，越近似
            while (Ta <= 2*Tj || Td <= 2*Tj)
                a_max = gamma*a_max;
                Tv = 0;
                Tj = a_max / j_max;
                Tj1 = Tj;
                Tj2 = Tj;
                delta = (a_max^4/j_max^2) + 2*(v_0^2 + v_1^2) + a_max*(4*(q_1 - q_0) - 2*(a_max/j_max)*(v_0 + v_1));
                Ta = ((power(a_max, 2)/j_max) - 2.0*v_0 + sqrt(delta)) / (2.0*a_max);
                Td = ((power(a_max, 2)/j_max) - 2.0*v_1 + sqrt(delta)) / (2.0*a_max);
                % 对Ta和Td进行讨论
                if Ta < 0 || Td < 0
                    if Ta < 0 && v_0 > v_1
                        % 只有减速段,Ta=0，结合公式（3.28）计算
                        Ta = 0;
                        Tj1 = 0;
                        Td = 2*(q_1 - q_0) / (v_0 + v_1);
                        Tj2 = (j_max*(q_1 - q_0) - sqrt(j_max*(j_max*power(q_1 - q_0, 2) + power(v_1 + v_0, 2)*(v_1 - v_0)))) / (j_max*(v_1 + v_0));
                        vlim = v_0;
                        calresult = [Ta, Tv, Td, Tj1, Tj2, q_0, q_1, v_0, v_1, vlim, j_max, j_min, sigma];
                        return
                    elseif Td < 0 && v_1 > v_0
                        % 只有加速段，Td=0，结合公式（3.29）计算
                        Td = 0;
                        Tj2 = 0;
                        Ta = 2*(q_1 - q_0) / (v_0 + v_1);
                        Tj1 = (j_max*(q_1 - q_0) - sqrt(j_max*(j_max*power(q_1 - q_0, 2)) - power(v_1 + v_0, 2)*(v_1 - v_0))) / (j_max*(v_1 + v_0));
                        a_lima = j_max * Tj1;
                        vlim = v_0 + (Ta - Tj1) * a_lima;
                        calresult = [Ta, Tv, Td, Tj1, Tj2, q_0, q_1, v_0, v_1, vlim, j_max, j_min, sigma];
                        return
                    end
                    
                elseif  Ta > 2*Tj && Td > 2*Tj
                    % 加速和减速阶段均能达到最大加速度
                    a_lima = j_max * Tj1;
                    vlim = v_0 + a_lima*(Ta - Tj);
                    calresult = [Ta, Tv, Td, Tj1, Tj2, q_0, q_1, v_0, v_1, vlim, j_max, j_min, sigma];
                    return
                end
            end
        end
    end
end
