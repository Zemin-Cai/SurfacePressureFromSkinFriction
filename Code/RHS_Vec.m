function [T] = RHS_Vec(Pressure_Field, A, B, C, D, E, F, h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construting the right hand size vector whose size is
% [(height-2)*(width-2),1]
%
% 04.18.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[height, width] = size(Pressure_Field);
T = zeros((height-2)*(width-2), 1);

for j = 1:1:(height-2)
    for i = 1:1:(width-2)
        if((i==1) && (j==1))
            T((i-1)*(height-2)+j) = 0.5*B(j,i)*Pressure_Field(j,i) + (A(j,i)*h^(-2) - 0.5*D(j,i)*h^(-1))*Pressure_Field(j+1,i) - 0.5*B(j,i)*Pressure_Field(j+2, i) + (C(j,i)*h^(-2) - 0.5*E(j,i)*h^(-1))*Pressure_Field(j, i+1) - 0.5*B(j,i)*Pressure_Field(j, i+2) - F(j,i);
        else
            if(i==(width-2) && j==(height-2))
                T((i-1)*(height-2)+j) = -0.5*B(j,i)*Pressure_Field(j,i+2) + (A(j,i)*h^(-2) + 0.5*D(j, i)*h^(-1))*Pressure_Field(j+1,i+2) + 0.5*B(j,i)*Pressure_Field(j+2, i+2) + (C(j,i)*h^(-2) + 0.5*E(j, i)*h^(-1))*Pressure_Field(j+2, i+1) - 0.5*B(j,i)*Pressure_Field(j+2, i) - F(j, i);   % ?
            else
                if(i==1 && j==(height-2))
                    T((i-1)*(height-2)+j) = 0.5*B(j,i)*Pressure_Field(j,i) + (A(j,i)*h^(-2) - 0.5*D(j,i)*h^(-1))*Pressure_Field(j+1,i) - 0.5*B(j,i)*Pressure_Field(j+2, i) + (C(j,i)*h^(-2)+0.5*E(j,i)*h^(-1))*Pressure_Field(j+2, i+1) + 0.5*B(j,i)*Pressure_Field(j+2, i+2) - F(j, i);
                else
                    if(i==(width-2) && j==1)
                        T((i-1)*(height-2)+j) = 0.5*B(j,i)*Pressure_Field(j,i) + (C(j,i)*h^(-2)-0.5*E(j,i)*h^(-1))*Pressure_Field(j,i+1) - 0.5*B(j,i)*Pressure_Field(j, i+2) + (A(j,i)*h^(-2)+0.5*D(j,i)*h^(-1))*Pressure_Field(j+1, i+2) + 0.5*B(j,i)*Pressure_Field(j+2, i+2) - F(j, i);
                    else
                        if(i==1 && j>1 && j<(height-2))
                            T((i-1)*(height-2)+j) = 0.5*B(j,i)*Pressure_Field(j,i) + (A(j,i)*h^(-2)-0.5*D(j,i)*h^(-1))*Pressure_Field(j+1,i) - 0.5*B(j,i)*Pressure_Field(j+2, i) - F(j, i);
                        else
                            if(i==(width-2) && j>1 && j<(height-2))
                                T((i-1)*(height-2)+j) = -0.5*B(j,i)*Pressure_Field(j,i+2) + (A(j,i)*h^(-2)+0.5*D(j,i)*h^(-1))*Pressure_Field(j+1,i+2) + 0.5*B(j,i)*Pressure_Field(j+2, i+2) - F(j, i);
                            else
                                if(j==1 && i>1 && i<(width-2))
                                    T((i-1)*(height-2)+j) = 0.5*B(j,i)*Pressure_Field(j,i) + (C(j,i)*h^(-2)-0.5*E(j,i)*h^(-1))*Pressure_Field(j,i+1) - 0.5*B(j,i)*Pressure_Field(j, i+2) - F(j, i);
                                else
                                    if(j==(height-2) && i>1 && i<(width-2))
                                        T((i-1)*(height-2)+j) = -0.5*B(j,i)*Pressure_Field(j+2,i) + (C(j,i)*h^(-2)+0.5*E(j,i)*h^(-1))*Pressure_Field(j+2,i+1) + 0.5*B(j,i)*Pressure_Field(j+2, i+2) - F(j,i);
                                    else
                                        T((i-1)*(height-2)+j) = - F(j,i);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end