% x = 2*binvar(4*NT,1)-1;
x = 2*binvar(2*NT,1)-1;
y = sdpvar(1,1);
% ����Լ��
% Constraints = [sum(x) <= 10, x(1)==0, 0.5 <= x(2) <= 1.5];	
% for i = 1:7
% 	Constraints = [Constraints, x(i) + x(i+1)  <= x(i+2) + x(i+3)];
% end
Constraints = [y<=GR*(x)/sqrt(2),y<=0];
options= sdpsettings;
options.solver = 'intlinprog';
Objective = -y;
% Objective = norm(y-H*x)^2;
%�������
sol = optimize(Constraints, Objective, options)
%���������־
if sol.problem == 0
	% ˵������ɹ�,��ʱչʾ�õ��Ľ��
% 	solution = P*round(value(x))
    solution = round(value(x))
else
	display('������');
	sol.info
	yalmiperror(sol.problem)
end

