function result = evaluateC(CF, ~, ~, ~, ~, ~)
A=diff(CF, 1, 1).^2;
result =sum(sum(A(1,:,:)));
end