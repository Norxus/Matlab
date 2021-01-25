% x = 2*binvar(4*NT,1)-1;
x = 2*binvar(2*NT,1)-1;
y = sdpvar(1,1);
% 定义约束
% Constraints = [sum(x) <= 10, x(1)==0, 0.5 <= x(2) <= 1.5];	
% for i = 1:7
% 	Constraints = [Constraints, x(i) + x(i+1)  <= x(i+2) + x(i+3)];
% end
Constraints = [y<=GR*(x)/sqrt(2),y<=0];
options= sdpsettings;
options.solver = 'intlinprog';
Objective = -y;
% Objective = norm(y-H*x)^2;
%求解问题
sol = optimize(Constraints, Objective, options)
%分析错误标志
if sol.problem == 0
	% 说明计算成功,此时展示得到的结果
% 	solution = P*round(value(x))
    solution = round(value(x))
else
	display('错了亲');
	sol.info
	yalmiperror(sol.problem)
end

