load("C:\Users\user\Downloads\iddata-03.mat");
plot(id),title("Datele de identificare")
figure,plot(val),title("Datele de validare")

identificare.u = id.u;
identificare.y = id.y;

validare.u = val.u;
validare.y = val.y;

na=3;
nb=3;
m =2;
ok_m1=1;
ok_m2=1;
ok_m3=1;
phi = [];
for k = 1:length(id.y)
   if(ok_m1 == 1)
    for index_coloana_grad1 = 1:na
        if(k - index_coloana_grad1<=0)
            phi(k,index_coloana_grad1) = 0;
        end
        if(k - index_coloana_grad1 > 0)
            phi(k,index_coloana_grad1) = identificare.y(k - index_coloana_grad1);
        end
    end
    for index_coloana_grad1 = na+1:na+nb
        if(k - (index_coloana_grad1 - na) <=0)
            phi(k,index_coloana_grad1) = 0;
        end
        if(k - (index_coloana_grad1 - na) > 0)
            phi(k,index_coloana_grad1) = identificare.u(k - (index_coloana_grad1 - na));
        end
    end
  end
    if(m == 1)
        index_coloana_grad1=index_coloana_grad1+1;
        phi(k,index_coloana_grad1)=1;
        ok_m2 = 0;
        ok_m3=0;
    end

   if(ok_m2==1 && ok_m1 == 1)
    index_coloana_grad2 =na+nb;
    for i = 1:na+nb
        for j = i:na+nb
             index_coloana_grad2 = index_coloana_grad2 + 1;
             phi(k,index_coloana_grad2) = phi(k,i) * phi(k,j);
        end
    end
   end
    
    if(m==2)
        index_coloana_grad2=index_coloana_grad2+1;
        phi(k,index_coloana_grad2)=1;
        ok_m3=0;
    end

  if(ok_m3 == 1 && ok_m2==1 && ok_m1==1)
    aux_vector = [];
    contor = 1;
    for index_aux = na+nb+1:index_coloana_grad2
        for index_aux2 = 1:na+nb
            aux_vector(contor) = phi(k,index_aux) * phi(k,index_aux2);
            contor = contor + 1;
        end
    end

    vec_nerepetat = [];
    indice_nerepetat = 0;
    lungime_aux = length(aux_vector);
    for t = 1:lungime_aux
        ok = 1;
        for t1 = t+1:lungime_aux
            if(aux_vector(t) == aux_vector(t1))
               ok = 0;
            end
        end
        if(ok == 1)
           indice_nerepetat = indice_nerepetat + 1;
           vec_nerepetat(indice_nerepetat) = aux_vector(t); 
        end
    end
  
    indice_nerepetat = 0;
    index_coloana_grad3 = index_coloana_grad2;
    for i = na+nb+1:index_coloana_grad2
        for j = 1:na+nb
            index_coloana_grad3 = index_coloana_grad3 + 1;
            indice_nerepetat = indice_nerepetat + 1;
            if(indice_nerepetat <= length(vec_nerepetat))
               phi(k,index_coloana_grad3) = vec_nerepetat(indice_nerepetat);
            end
        end
    end
  end
    if(m==3)
        index_coloana_grad3=index_coloana_grad3+1;
        phi(k,index_coloana_grad3)=1;
    end
end

theta = phi\identificare.y;
yH_id = phi * theta;
yH_id(1) = identificare.y(1);
figure,plot(yH_id,"b"),title("Pt datele de Identificare-Predictie")
hold on;
plot(identificare.y,"r")

MSE_identificare_predictie = 0;
MSE_identificare_predictie = sum((yH_id - identificare.y).^2)/length(id.y)

ok_m1_val=1;
ok_m2_val=1;
ok_m3_val=1;
phi_validare = [];
for k = 1:length(val.y)
   if(ok_m1_val==1)
    for index_coloana_grad1 = 1:na
        if(k - index_coloana_grad1<=0)
            phi_validare(k,index_coloana_grad1) = 0;
        end
        if(k - index_coloana_grad1 > 0)
            phi_validare(k,index_coloana_grad1) = validare.y(k - index_coloana_grad1);
        end
    end
    for index_coloana_grad1 = na+1:na+nb
        if(k - (index_coloana_grad1 - na) <=0)
            phi_validare(k,index_coloana_grad1) = 0;
        end
        if(k - (index_coloana_grad1 - na) > 0)
            phi_validare(k,index_coloana_grad1) = validare.u(k - (index_coloana_grad1 - na));
        end
    end
   end

   if(m==1)
       index_coloana_grad1=index_coloana_grad1+1;
       phi_validare(k,index_coloana_grad1)=1;
       ok_m2_val=0;
       ok_m3_val=0;
   end
   if(ok_m2_val==1 && ok_m1_val==1)
    index_coloana_grad2 = na + nb;
    for i = 1:na+nb
        for j = i:na+nb
             index_coloana_grad2 = index_coloana_grad2 + 1;
             phi_validare(k,index_coloana_grad2) = phi_validare(k,i) * phi_validare(k,j);
        end
    end
   end
    if(m==2)
        index_coloana_grad2=index_coloana_grad2+1;
        phi_validare(k,index_coloana_grad2)=1;
        ok_m3_val=0;
    end

   if(ok_m3_val == 1 && ok_m2_val==1 && ok_m3_val==1)
    aux_vector_val = [];
    contor = 1;
    for index_aux = na+nb+1:index_coloana_grad2
        for index_aux2 = 1:na+nb
            aux_vector_val(contor) = phi_validare(k,index_aux) * phi_validare(k,index_aux2);
            contor = contor + 1;
        end
    end

    vec_nerepetat = [];
    indice_nerepetat = 0;
    lungime_aux = length(aux_vector_val);
    for t = 1:lungime_aux 
        ok = 1;
        for t1 = t+1:lungime_aux
            if(aux_vector_val(t) == aux_vector_val(t1))
               ok = 0;
            end
        end
        if(ok == 1)
           indice_nerepetat = indice_nerepetat + 1;
           vec_nerepetat(indice_nerepetat) = aux_vector_val(t); 
        end
    end

    indice_nerepetat = 0;
    index_coloana_grad3 = index_coloana_grad2;
    for i = na+nb+1:index_coloana_grad2
        for j = 1:na+nb
            index_coloana_grad3 = index_coloana_grad3 + 1;
            indice_nerepetat = indice_nerepetat + 1;
             if(indice_nerepetat <= length(vec_nerepetat))
                phi_validare(k,index_coloana_grad3) = vec_nerepetat(indice_nerepetat);
             end
        end
    end
   end
   if(m==3)
       index_coloana_grad3=index_coloana_grad3+1;
       phi_validare(k,index_coloana_grad3)=1;
   end
end

figure;

yH_val = phi_validare * theta;
yH_val(1) = validare.y(1);
plot(yH_val,"b"),title("Pt datele de Validare-Predictie")
hold on;
plot(validare.y,"r")

MSE_validare_predictie = 0;
MSE_validare_predictie = sum((yH_val - validare.y).^2)/length(val.y)

%simularea modelului pe datele de identificare
y_identificare_simulare=[];
dimensiune=length(yH_id);
y_identificare_simulare(1)=identificare.y(1);
for index_linie=2:dimensiune
   if(m>=1)
    y_identificare_simulare(index_linie)=0;
    for coloana=1:na
        if((index_linie-coloana)>0)
            y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+y_identificare_simulare(index_linie-coloana)*theta(coloana);
        else
            y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+0;
        end
    end
    for coloana_2=(na+1):(na+nb)
        if(index_linie-(coloana_2-na)>0)
            y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+identificare.u(index_linie-coloana_2+na)*theta(coloana_2);
        else
            y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+0;
        end
    end
  end
    if(m==1)
       y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+theta(coloana_2+1);
    end
  if(m>1)
    indice_theta=coloana_2+1;
    for indice_linie=1:na+nb
        for indice_coloana=indice_linie:na+nb
            if(indice_coloana<=na && indice_linie<=na)
                aux=0;
                if((index_linie-indice_linie)>=1 && (index_linie-indice_coloana)>=1)
                    aux=y_identificare_simulare(index_linie-indice_linie)*y_identificare_simulare(index_linie-indice_coloana);
                    y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+aux*theta(indice_theta);
                    indice_theta=indice_theta+1;
                else
                    y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+0*theta(indice_theta);
                    indice_theta=indice_theta+1;
                end
            elseif(indice_coloana>na && indice_linie<=na)
                aux1=0;
                if((index_linie-indice_linie)>=1 && (index_linie-(indice_coloana-na))>=1)
                    aux1=y_identificare_simulare(index_linie-indice_linie)*identificare.u(index_linie-(indice_coloana-na));
                    y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+aux1*theta(indice_theta);
                    indice_theta=indice_theta+1;
                else
                    y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+0*theta(indice_theta);
                    indice_theta=indice_theta+1;
                end
            else
                aux2=0;
                if((index_linie-(indice_linie-na))>=1 && (index_linie-(indice_coloana-na))>=1)
                    aux2=identificare.u(index_linie-(indice_linie-na))*identificare.u(index_linie-(indice_coloana-na));
                    y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+aux2*theta(indice_theta);
                    indice_theta=indice_theta+1;
                else
                    y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+0*theta(indice_theta);
                    indice_theta=indice_theta+1;
                end
            end
       end
    end
     y_identificare_simulare(index_linie)=y_identificare_simulare(index_linie)+theta(indice_theta);
  end
end
MSE_identificare_simulare=0;
MSE_identificare_simulare=sum((y_identificare_simulare'-identificare.y).^2)/length(identificare.y)
figure,plot(y_identificare_simulare,'b'),title('Datele de identificare-Simulare')
hold on,plot(identificare.y,'r')
%simularea modelului pe datele de validare

y_validare_simulare=[];
dimensiune1=length(yH_val);
y_validare_simulare(1)=validare.y(1);
for index_linie=2:dimensiune1
  if(m>=1)
    y_validare_simulare(index_linie)=0;
    for coloana=1:na
        if((index_linie-coloana)>0)
            y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+y_validare_simulare(index_linie-coloana)*theta(coloana);
        else
            y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+0;
        end
    end
    for coloana_2=(na+1):(na+nb)
        if(index_linie-(coloana_2-na)>0)
            y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+validare.u(index_linie-coloana_2+na)*theta(coloana_2);
        else
            y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+0;
        end
    end
  end
   if(m==1)
       y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+theta(coloana_2+1);
   end
if(m>1)
    indice_theta1=coloana_2+1;
    for indice_linie=1:na+nb
        for indice_coloana=indice_linie:na+nb
            if(indice_coloana<=na && indice_linie<=na)
                auxiliar=0;
                if((index_linie-indice_linie)>=1 && (index_linie-indice_coloana)>=1)
                    auxiliar=y_validare_simulare(index_linie-indice_linie)*y_validare_simulare(index_linie-indice_coloana);
                    y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+auxiliar*theta(indice_theta1);
                    indice_theta1=indice_theta1+1;
                else
                    y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+0;
                    indice_theta1=indice_theta1+1;
                end
            elseif(indice_coloana>na && indice_linie<=na)
                auxiliar1=0;
                if((index_linie-indice_linie)>=1 && (index_linie-(indice_coloana-na))>=1)
                    auxiliar1=y_validare_simulare(index_linie-indice_linie)*validare.u(index_linie-(indice_coloana-na));
                    y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+auxiliar1*theta(indice_theta1);
                    indice_theta1=indice_theta1+1;
                else
                    y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+0;
                    indice_theta1=indice_theta1+1;
                end
            else
                auxiliar2=0;
                if((index_linie-(indice_linie-na))>=1 && (index_linie-(indice_coloana-na))>=1)
                    auxiliar2=validare.u(index_linie-(indice_linie-na))*validare.u(index_linie-(indice_coloana-na));
                    y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+auxiliar2*theta(indice_theta1);
                    indice_theta1=indice_theta1+1;
                else
                    y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+0;
                    indice_theta1=indice_theta1+1;
                end
            end
       end
    end
    y_validare_simulare(index_linie)=y_validare_simulare(index_linie)+theta(indice_theta1);
 end
end

MSE_validare_simulare=sum((y_validare_simulare'-validare.y).^2)/length(validare.y)
figure,plot(y_validare_simulare,'b'),title('Datele de validare-Simulare')
hold on,plot(validare.y,'r')
