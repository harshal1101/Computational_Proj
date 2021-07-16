clc
clear all
Tc = 324.7;  %critical temp HCL in kelvin
Pc = 83.10; %critical pressure HCL in bar
R = 0.08314; % R in L bar (K-1) (mol-1)
a = (0.42748*(R*R*(Tc.^(2.5))))/Pc;
b = (0.08664*R*Tc)/Pc;   % b = (0.08664*RTc/Pc)
p_a = 4.57389; %Antoine A
p_b = 868.358; %Antoine B
p_c = 1.754;   %Antoine C
T = [298,300,302,304,306,308,310,312,314,316,318,320,322,324,324.7,326];
tolerance = 0.001;
%i=7;
v = linspace(0.03,3,1000);
for i=1:length(T)
c1 = (R*T(i));
c2 = (b*R*T(i) - (a/sqrt(T(i))));
c3 = (a*b)/sqrt(T(i));
d1  = b*b;
p = @(v) ((c1.*v.*v)+(c2.*v)+c3)./((v.*v.*v) -(d1.*v));
f = p(v);
p(2);
xlabel('V in L/mol','FontSize',40);
ylabel('P in bar','FontSize',40);
xlim([0 0.4])
ylim([20 100])
plot(v,f);
%axis([0.00002 0.0006 x y]);
hold on;
%T  = T+2;
end
function [n_roots,rho_l,rho_v,diff] = sol(psat,R,a,b,T,i)
  f1 = @(y) (R*T(i))/y;
  f2 = @(y) ((b*R*T(i)/y)-(a/(sqrt(T(i))*y))+(b*b));
  f3 = @(y) (a*b/(y*sqrt(T(i))));
  x = f1(psat);
  y = f2(psat);
  z=  f3(psat);
  poly = [1 -x -y -z];
  ro = roots(poly);
  ro = ro(imag(ro)==0);
  n_roots = length(ro);
  rho_l = (1./min(ro));
  rho_v = (1./max(ro));
  t2 = (a/(b*R*(T(i).^(3/2))));
  mu = @(x) (-log(1-x)-(t2*log(1+x)) + (x/(1-x)) -(t2*x/(1+x)) - log(psat*b/(x*R*T(i))));
  rho_lf = b*rho_l;
  rho_vf = b*rho_v;
  lhs =  mu(rho_lf);
  rhs = mu(rho_vf);
  diff = lhs -rhs;
endfunction
fpsat = [];
fvl =[];
fvv = [];
T = [300,302,304,306,308,310,312,314,316,318,320,322,324,324.7];
%T=324.7
for j=1:length(T)
  
  if j<=7
    calc= p_a -  (p_b)/(T(j)+p_c);
    psat = 10.^(calc)
  elseif j<length(T)-2
    psat = fpsat(j-1)+2
  elseif j<length(T)
    psat = fpsat(j-1)+2.5
  elseif j==length(T)
    psat = fpsat(j-1)+1
  end
  
  mini=100;
  %psat = 83.126;
  [nroots,rho_l,rho_v,diff] = sol(psat,R,a,b,T,j)
  disp(" ");
  cnt = 1;
  if nroots==3   %root of psat
    if abs(diff)<tolerance
      vl = 1./rho_l;
      vv = 1./rho_v;
      %plot([vl,vv],[psat,psat],'-o')
      fpsat(j) =psat;
      fvl(j) =vl;
      fvv(j) =vv;
      hold on;
    else   %if diff>tolerance
      %if cnt==1  %for first iteration taking two values (psat-10) and (psat+10) respectively
      [nrootsl,rho_ll,rho_vl,diffl] = sol(psat-10,R,a,b,T,j); %psat-10
      disp(" ");
      [nrootsu,rho_lu,rho_vu,diffu] = sol(psat+10,R,a,b,T,j); %psat+10      
      psatl = psat-10;psatu= psat+10;
      if nrootsl==1&&nrootsu==1
        psatl = (2*psat-10)/2; %lower guess value of pressure
        psatu = (2*psat+10)/2; %higher guess value of pressure
        [nrootsl,rho_ll,rho_vl,diffl] = sol(psatl,R,a,b,T,j);
        [nrootsu,rho_lu,rho_vu,diffu] = sol(psatu,R,a,b,T,j);
      elseif nrootsl==1
        %[nrootsl,rho_ll,rho_vl,diffl] = [nroots,rho_l,rho_v,diff]; %lower guess value
        nrootsu = nroots;
        rho_lu =rho_l;
        rho_vu = rho_v;
        diffu = diff;        
        psatl =psat;
      elseif nrootsu==1
        %[nroots,rho_ll,rho_vl,diffl] = [nroots,rho_l,rho_v,diff];
        nrootsu = nroots;
        rho_lu =rho_l;
        rho_vu = rho_v;
        diffu = diff;
        psatu = psat;
        fprintf("diffu diff %f %f\n",diffu,diff);
      else
        psatr = (psatu) - ((diffu)*(psatl-psatu)/(diffl-diffu)); %false-position
        [nrootsr,rho_lr,rho_vr,diffr] = sol(psatr,R,a,b,T,j)        
        if (diffl*diffr)<0
          %[nrootsu,rho_lu,rho_vu,diffu] = [n_rootsr,rho_lr,rho_vr,diffr];
          nrootsu = nrootsr;
          rho_lu =rho_lr;
          rho_vu = rho_vr;
          diffu = diffr;
          psatu =psatr; 
          fprintf(" i m in diffl and diffr \n");
        else
          %[nrootsl,rho_ll,rho_vl,diffl] = [n_rootsr,rho_lr,rho_vr,diffr];
          nrootsl = nrootsr;
          rho_ll =rho_lr;
          rho_vl = rho_vr;
          diffl = diffr;
          psatl =psatr; 
          psatl =psatr; 
          fprintf(" i m in diffu and diffr \n");
        end %if diffl*diffr<0
      end %if nroots==1&&nrootsu==1
      flag=0;
      for k=1:100
        %{
        if abs(diffu)<tolerance
          vlu = 1./rho_lu
          vvu = 1./rho_vu
          %plot([vlu,vvu],[psatu,psatu],'-o');
          hold on;
          fprintf(" diffu vlu vvu psatu %f %f %f %f\n",diffu,vlu,vvu,psatu);
          fpsat(j) = psatu;
          fvl(j) =vlu;
          fvv(j) =vvu;
          break;
        elseif abs(diffl)<tolerance
          vll = 1./rho_ll
          vvl = 1./rho_vl
          %plot([vll,vvl],[psatl,psatl],'-o');
          hold on;
          fprintf(" diffl vll vvl psatl %f %f %f %f \n",diffl,vll,vvl,psatl);
          fpsat(j) = psatl;
          fvl(j) =vll;
          fvv(j) =vvl;
          break;
        end 
        %}
        if nrootsl==1&&nrootsu==1
          psatl = (3*psatl+psatu)/4; %lower guess value of pressure
          psatu = (psatl+(3*psatu))/4 ; %higher guess value of pressure
          [nrootsl,rho_ll,rho_vl,diffl] = sol(psatl,R,a,b,T,j);
          [nrootsu,rho_lu,rho_vu,diffu] = sol(psatu,R,a,b,T,j);
        elseif nrootsl==1
          %[nrootsl,rho_ll,rho_vl,diffl] = [nroots,rho_l,rho_v,diff]; %lower guess value
          nrootsu = nroots;
          rho_lu =rho_l;
          rho_vu = rho_v;
          diffu = diff;        
          psatl =psat;
        elseif nrootsu==1
          %[nroots,rho_ll,rho_vl,diffl] = [nroots,rho_l,rho_v,diff];
          nrootsu = nroots;
          rho_lu =rho_l;
          rho_vu = rho_v;
          diffu = diff;
          psatu = psat;
        else
          psatr = (psatu) - ((diffu)*(psatl-psatu)/(diffl-diffu)); %false-position
          [nrootsr,rho_lr,rho_vr,diffr] = sol(psatr,R,a,b,T,j);
          if (diffl*diffr)<0
            %[nrootsu,rho_lu,rho_vu,diffu] = [n_rootsr,rho_lr,rho_vr,diffr];
            nrootsu = nrootsr;
            rho_lu =rho_lr;
            rho_vu = rho_vr;
            diffu = diffr;
            psatu =psatr; 
          else
            %[nrootsl,rho_ll,rho_vl,diffl] = [n_rootsr,rho_lr,rho_vr,diffr];
            nrootsl = nrootsr;
            rho_ll =rho_lr;
            rho_vl = rho_vr;
            diffl = diffr;
            psatl =psatr; 
            psatl =psatr; 
          end %if diffl*diffr<0
        end   %nroots==1  
        if abs(diffu)<tolerance&&nrootsu==3
          vlu = 1./rho_lu;
          vvu = 1./rho_vu;
          %plot([vlu,vvu],[psatu,psatu],'-o');
          hold on;
          fprintf(" diffu vlu vvu psatu %f %f %f %f\n",diffu,vlu,vvu,psatu);
          fpsat(j) = psatu;
          fvl(j) =vlu;
          fvv(j) =vvu;
          flag=1;
          break;
        elseif abs(diffl)<tolerance&&nrootsl==3
          vll = 1./rho_ll
          vvl = 1./rho_vl
          %plot([vll,vvl],[psatl,psatl],'-o');
          hold on;
          fprintf(" diffl vll vvl psatl %f %f %f %f \n",diffl,vll,vvl,psatl);
          fpsat(j) = psatl;
          fvl(j) =vll;
          fvv(j) =vvl;
          flag=1;
          break;
        end  
        %fprintf("diffl and diffu are %f %f\n",diffl,diffu);
        %fprintf("psat, vl and vg are %f %f %f\n",psatu,rho_lu,rho_vu);
        disp(k);
        
        if abs(diffl)<abs(diffu)
          if diffl==0
            if abs(diffu)<mini
              mini = diffu;
              minivl = 1./rho_lu;
              minivv = 1./rho_vu;
              minp = psatu;
            end %abs(diffl)<mini
          else
            if abs(diffl)<mini
              mini = diffl;
              minivl = 1./rho_ll;
              minivv = 1./rho_vl;
              minp = psatl;
            end %abs(diffu)<mini
          end %diffu==0 
        end
        if abs(diffu)<abs(diffl)
          if diffu==0
            if abs(diffl)<mini
              mini = diffl;
              minivl = 1./rho_ll;
              minivv = 1./rho_vl;
              minp = psatl;
            end %abs(diffl)<mini
          else
            if abs(diffu)<mini
              mini = diffu;
              minivl = 1./rho_lu;
              minivv = 1./rho_vu;
              minp = psatu;
            end %abs(diffu)<mini
          end %diffu==0
        end %abs(diffu)<abs(diffl)
        
      end % for k=1:100
      if flag==0
        fpsat(j) = minp;
        fvl(j) = minivl;
        fvv(j) =minivv;
      end      
    end % end of diff<tolerance
  else
    [nrootsl,rho_ll,rho_vl,diffl] = sol(psat-0.5,R,a,b,T,j); %psat-10
    disp(" ");
    [nrootsu,rho_lu,rho_vu,diffu] = sol(psat+0.5,R,a,b,T,j); %psat+10      
    psatl = psat-0.5;psatu= psat+0.5;
    if nrootsl==1&&nrootsu==1
      psatl = (2*psat-0.5)/2 %lower guess value of pressure
      psatu = (2*psat+0.5)/2 %higher guess value of pressure
      [nrootsl,rho_ll,rho_vl,diffl] = sol(psatl,R,a,b,T,j);
      [nrootsu,rho_lu,rho_vu,diffu] = sol(psatu,R,a,b,T,j);
      fpsat(j) =psatu;
      fvl(j) = 1./rho_l;
      fvv(j) = 1./rho_v;
    elseif nrootsl==1
      %[nrootsl,rho_ll,rho_vl,diffl] = [nroots,rho_l,rho_v,diff]; %lower guess value
     % nrootsu = nroots;
      %rho_ll =rho_l;
      %rho_vl = rho_v;
      %diffl = diff;        
      %psatl =psat;
      fpsat(j) =psatu;
      fvl(j) = 1./rho_lu;
      fvv(j) = 1./rho_vu;
     % printf("diffu diff %f %f\n",diffl,diff);
    elseif nrootsu==1
      %[nroots,rho_ll,rho_vl,diffl] = [nroots,rho_l,rho_v,diff];
      %nrootsu = nroots;
      %rho_lu =rho_l;
      %rho_vu = rho_v;
      %diffu = diff;
      %psatu = psat;
      fpsat(j) =psatl;
      fvl(j) = 1./rho_ll;
      fvv(j) = 1./rho_vl;
      %printf("diffu diff %f %f\n",diffu,diff);
   
    end %nroots=1
  end %end of nroots=3
  
end
disp(fpsat);
disp(fvl);
disp(fvv);
pl = 1:length(T)
%plot(fvl(pl),fpsat(pl))
%plot(fvv(pl),fpsat(pl))
plz = polyfit(fvl,fpsat,6);
y1 = polyval(plz,fvl);
plot(fvl,y1,'LineWidth',8);
hold on;
plz2 = polyfit(fvv,fpsat,6);
y2 = polyval(plz2,fvv);
plot(fvv,y2,'LineWidth',8);







