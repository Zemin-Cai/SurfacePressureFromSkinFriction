function [M] = CoeffMatrix(height, width, Para_1, Para_2, Para_3, Para_4, Para_5, h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construting the coefficient matrix with size 
%  [(height-2)*(width-2),(height-2)*(width-2)
%
% 04.17.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = zeros((height-2)*(width-2), (height-2)*(width-2));
for j = 1:1:(height-2)
    for i = 1:1:(width-2)
        if((i==1) && (j==1))
            M(i, (i-1)*(height-2)+1) = (2*Para_1(j,i)*h^(-2) + 2*Para_3(j,i)*h^(-2));
            M(i, (i-1)*(height-2)+2) = -Para_3(j,i)*h^(-2) - 0.5 * Para_5(j, i)*h^(-1);
            M(i, i*(height-2)+1) = -Para_1(j,i)*h^(-2) - 0.5*Para_4(j, i)*h^(-1);
            M(i, i*(height-2)+2) = -0.5*Para_2(j,i)*h^(-2);
        else
            if(i==(width-2) && j==(height-2))
                M((i-1)*(height-2)+j,(i-2)*(height-2)+j-1) = -0.5*Para_2(j,i);
                M((i-1)*(height-2)+j,(i-2)*(height-2)+j) = -Para_1(j,i)*h^(-2) + 0.5*Para_4(j, i)*h^(-1);
                M((i-1)*(height-2)+j,(i-1)*(height-2)+j-1) = -Para_3(j,i)*h^(-2) + 0.5*Para_5(j, i)*h^(-1);
                M((i-1)*(height-2)+j,(i-1)*(height-2)+j) = 2*Para_1(j,i)*h^(-2) + 2*Para_3(j,i)*h^(-2);
            else
                if(i==1 && j==(height-2))
                    M((i-1)*(height-2)+j,(i-1)*(height-2)+j-1) = -Para_3(j,i)*h^(-2) + 0.5*Para_5(j, i)*h^(-1);
                    M((i-1)*(height-2)+j,(i-1)*(height-2)+j) = 2*Para_1(j,i)*h^(-2) + 2*Para_3(j,i)*h^(-2);
                    M((i-1)*(height-2)+j,i*(height-2)+j-1) = 0.5*Para_2(j,i);
                    M((i-1)*(height-2)+j,i*(height-2)+j) = -Para_1(j,i)*h^(-2) - 0.5*Para_4(j,i)*h^(-1);
                else
                    if(i==(width-2) && j==1)
                        M((i-1)*(height-2)+j,(i-2)*(height-2)+j) = -Para_1(j,i)*h^(-2) + 0.5*Para_4(j, i)*h^(-1);
                        M((i-1)*(height-2)+j,(i-2)*(height-2)+j+1) = 0.5*Para_2(j,i);
                        M((i-1)*(height-2)+j,(i-1)*(height-2)+j) = 2*Para_1(j,i)*h^(-2) + 2*Para_3(j,i)*h^(-2);
                        M((i-1)*(height-2)+j,(i-1)*(height-2)+j+1) = -Para_3(j,i)*h^(-2) - 0.5*Para_5(j, i)*h^(-1);
                    else
                        if(i==1 && j>1 && j<(height-2))
                            M((i-1)*(height-2)+j,(i-1)*(height-2)+j-1) = -Para_3(j,i)*h^(-2) + 0.5*Para_5(j, i)*h^(-1);
                            M((i-1)*(height-2)+j,(i-1)*(height-2)+j) = 2*Para_1(j,i)*h^(-2) + 2*Para_3(j,i)*h^(-2);
                            M((i-1)*(height-2)+j,(i-1)*(height-2)+j+1) = -Para_3(j,i)*h^(-2) - 0.5*Para_5(j, i)*h^(-1);
                            M((i-1)*(height-2)+j,i*(height-2)+j-1) = 0.5*Para_2(j,i);
                            M((i-1)*(height-2)+j,i*(height-2)+j) = -Para_1(j,i)*h^(-2) - 0.5*Para_4(j, i)*h^(-1);
                            M((i-1)*(height-2)+j,i*(height-2)+j+1) = -0.5*Para_2(j,i);
                        else
                            if(i==(width-2) && j>1 && j<(height-2))
                                M((i-1)*(height-2)+j,(i-2)*(height-2)+j-1) = -0.5*Para_2(j,i);
                                M((i-1)*(height-2)+j,(i-2)*(height-2)+j) = -Para_1(j,i)*h^(-2) + 0.5*Para_4(j,i)*h^(-1);
                                M((i-1)*(height-2)+j,(i-2)*(height-2)+j+1) = 0.5*Para_2(j,i);
                                M((i-1)*(height-2)+j,(i-1)*(height-2)+j-1) = -Para_3(j,i)*h^(-2) + 0.5*Para_5(j, i)*h^(-1);
                                M((i-1)*(height-2)+j,(i-1)*(height-2)+j) = 2*Para_1(j,i)*h^(-2) + 2*Para_3(j,i)*h^(-2);
                                M((i-1)*(height-2)+j,(i-1)*(height-2)+j+1) = -Para_3(j,i)*h^(-2) - 0.5*Para_5(j, i)*h^(-1);
                            else
                                if(j==1 && i>1 && i<(width-2))
                                    M((i-1)*(height-2)+j,(i-2)*(height-2)+j) = -Para_1(j,i)*h^(-2) + 0.5*Para_4(j, i)*h^(-1);
                                    M((i-1)*(height-2)+j,(i-2)*(height-2)+j+1) = 0.5*Para_2(j,i);
                                    M((i-1)*(height-2)+j,(i-1)*(height-2)+j) = 2*Para_1(j,i)*h^(-2) + 2*Para_3(j,i)*h^(-2);
                                    M((i-1)*(height-2)+j,(i-1)*(height-2)+j+1) = -Para_3(j,i)*h^(-2) - 0.5*Para_5(j, i)*h^(-1);
                                    M((i-1)*(height-2)+j,i*(height-2)+j) = -Para_1(j,i)*h^(-2) - 0.5*Para_4(j, i)*h^(-1);
                                    M((i-1)*(height-2)+j,i*(height-2)+j+1) = -0.5*Para_2(j,i);
                                else
                                    if(j==(height-2) && i>1 && i<(width-2))
                                        M((i-1)*(height-2)+j,(i-2)*(height-2)+j-1) = -0.5*Para_2(j,i);
                                        M((i-1)*(height-2)+j,(i-2)*(height-2)+j) = -Para_1(j,i)*h^(-2) + 0.5*Para_4(j, i)*h^(-1);
                                        M((i-1)*(height-2)+j,(i-1)*(height-2)+j-1) = -Para_3(j,i)*h^(-2) + 0.5*Para_5(j, i)*h^(-1);
                                        M((i-1)*(height-2)+j,(i-1)*(height-2)+j) = 2*Para_1(j,i)*h^(-2) + 2*Para_3(j,i)*h^(-2);
                                        M((i-1)*(height-2)+j,(i)*(height-2)+j-1) = 0.5*Para_2(j,i);
                                        M((i-1)*(height-2)+j,(i)*(height-2)+j) = -Para_1(j,i)*h^(-2) - 0.5*Para_4(j, i)*h^(-1);
                                    else
                                        M((i-1)*(height-2)+j,(i-2)*(height-2)+j-1) = -0.5*Para_2(j,i);
                                        M((i-1)*(height-2)+j,(i-2)*(height-2)+j) = -Para_1(j,i)*h^(-2) + 0.5*Para_4(j, i)*h^(-1);
                                        M((i-1)*(height-2)+j,(i-2)*(height-2)+j+1) = 0.5*Para_2(j,i);
                                        M((i-1)*(height-2)+j,(i-1)*(height-2)+j-1) = -Para_3(j,i)*h^(-2) + 0.5*Para_5(j, i)*h^(-1);
                                        M((i-1)*(height-2)+j,(i-1)*(height-2)+j) = 2*Para_1(j,i)*h^(-2) + 2*Para_3(j,i)*h^(-2);
                                        M((i-1)*(height-2)+j,(i-1)*(height-2)+j+1) = -Para_3(j,i)*h^(-2) - 0.5*Para_5(j, i)*h^(-1);
                                        M((i-1)*(height-2)+j,(i)*(height-2)+j-1) = 0.5*Para_2(j,i);
                                        M((i-1)*(height-2)+j,(i)*(height-2)+j) = -Para_1(j,i)*h^(-2) - 0.5*Para_4(j, i)*h^(-1);
                                        M((i-1)*(height-2)+j,(i)*(height-2)+j+1) = -0.5*Para_2(j,i);
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