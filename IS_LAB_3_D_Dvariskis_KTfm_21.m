%Domantas Dvariškis
%KTfm-21
%3 Laboratorinis darbas pagrindinė užd.
%Netiesinio aproksimatoriaus,spindulio tipo bazinių funkcijų tinklu grįsto mokymo algoritmas

clear all
clc

%1. Duomenų paruošimas mokymui
x = [0.1:1/22:1];
d = (1 + 0.6*sin(2*pi*x/0.7) + 0.3*sin(2*pi*x))./2; %Išėjime norimas atsakas
figure(1);
plot(x,d,'b*')
grid on;

% pradinių ryšių svorių generavimas
w1 = rand(1)*0.1;
w2 = rand(1)*0.1;
w0 = rand(1)*0.1;
c1 = 0.19; 
r1 = 0.2;
c2 = 0.91;
r2 = 0.21; 

%F = exp(-(x-c)^2/(2*r^2));

n = 0.15; %Mokymo žingsnis
% Mokymas
for i = 1:100 %Tinklo apmokymo pakartojimai
    for i = 1:length(x)
        %Paslėptojo sluoksnio lokalios išeities skaičiavimas
        F1 = exp(-(x(i)-c1)^2/(2*r1^2));
        F2 = exp(-(x(i)-c2)^2/(2*r2^2));
        % aktyvavimo funkcijos pritaikymas 
        y = F1*w1 + F2*w2 + w0;
        % palyginama su norimu atsaku ir skaičiavimo klaida
        e = d(i) - y;
        % tinklo ryšių svorių atnaujinimas
        w1 = w1 + n*e*F1;
        w2 = w2 + n*e*F2;
        w0 = w0 + n*e;       
    end
end
    
%Gaunamos tinklo reikšmės po apmokymo
Y = zeros(1,length(d));
for i = 1:length(x)
    %Paslėptojo sluoksnio lokalios išeities skaičiavimas
    F1 = exp(-(x(i)-c1)^2/(2*r1^2));
    F2 = exp(-(x(i)-c2)^2/(2*r2^2));
    % aktyvavimo funkcijos pritaikymas 
    y = F1*w1 + F2*w2 + w0;
    %Saugome tarpines reikšmes
    Y(i) = y;
end

hold on 
plot(x,Y,'r+')
hold off;
