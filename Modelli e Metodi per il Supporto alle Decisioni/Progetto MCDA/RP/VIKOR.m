%--------------------------------------------------------------
%--------------------------------------------------------------
% METODO VIKOR
%--------------------------------------------------------------
%--------------------------------------------------------------

clear
% LETTURA DATI DI INPUT
    X = dlmread('Mlocal.txt');
    disp('Dimensions: (m alternatives, n criteria)')
    [m,n] = size(X)
    disp('Criteria directions (1: max, -1: min)')
    Dir = dlmread('Dir.txt')
    disp('Criteria weights ')
    W = dlmread('W.txt')
    v=0.5
    delta =1/(m-1);
% NORMALIZZAZIONE [0,1] con TRASFORMAZIONE CRITERI (min -> max).
% MN = Normalized Matrix.
    for j = 1:n
        if (Dir(j)==1)
            val_min = min(X(:,j));
            val_max = max(X(:,j));
        else
             val_max = min(X(:,j));
             val_min = max(X(:,j)); 
        end
        for i = 1:m 
            MN(i,j) = (X(i,j)-val_min)/(val_max-val_min);
        end 
    end
% MATRICE PESATA V (V(i,j)=w(j) * MN(i,j))
    for j=1:n
      V(:,j) = W(j) * MN(:,j);
    end
% CALCOLO SOLUZIONE IDEALE.
    I=W; 
% CALCOLO DELLE DISTANZE S ED R.
% S: distanza norma_1 dall'ideale.
% R: distanza norma_inf dall'ideale.
    for i=1:m
      S(i) = norm(I-V(i,:),1); 
      R(i) = norm(I-V(i,:),inf);
    end
% CALCOLO DISTANZE QS E QR NORMALIZZATE IN [0,1]
    Smin = min(S);
    Smax = max(S);
    for i = 1:m 
        QS(i) = (S(i)-Smin)/(Smax-Smin);
    end
    Rmin = min(R);
    Rmax = max(R);
    for i = 1:m 
         QR(i) = (R(i)-Rmin)/(Rmax-Rmin);
    end 
% CALCOLO Q COME AGGREGAZIONE QR e QS
    for i = 1:m 
         Q(i) = (v*QS(i))+((1-v)*QR(i));
    end     
% ORDINAMENTO INDICE Q e STAMPA CLASSIFICA
    disp('VIKOR Ranking (increasing index Q)')
    Ranking = sortrows( [ Q' (1:m)' ] , 'ascend' );
    for i = 1:m
        st = sprintf("%2d) \t%2d \t%f", i, Ranking(i,2), Ranking(i,1) );
        disp(st)
    end
    disp(' ')
% CALCOLO A' ed A''
% A' -> vincitore
    Q_primo = Ranking(1,1); %VALORE DI A'
    indQ_primo = Ranking(1,2); %INDICE DI A'
% A''-> secondo vincitore
    Q_secondo = Ranking(2,1); %VALORE DI A''
    indQ_secondo = Ranking(2,2); %INDICE DI A''
% VERIFICO C1 e C2
% C1 -> Q(a'')-Q(a') >= delta
    C1=false;
    if(Q_secondo-Q_primo>=delta)
        C1=true;
    end
% C2 -> A' vincitore secondo classifica QS o QR
C2=false;
    if(QR(indQ_primo)==min(QR) || QS(indQ_primo)==min(QS))
        C2=true;
    end
    
% COMPROMISE SOLUTION
% Verificate C1 e C2 -> a' vincitore
    if(C1 && C2)
        CS = indQ_primo;
        disp('Sono verificate sia le condizioni C1 che C2')
        disp('CS (compromise solution)= ');
        disp(CS);
    end
% Verificata C1 ma non C2 -> {a',a''}
    if(~C2 && C1)
        disp('Verificata la condizione C1, ma non la condizione C2')
        CS = [indQ_primo;indQ_secondo];
        disp('CS (compromise solution)= ');
        disp (CS)
    end
% Non verificata C1 -> {a(i):Q(i)-Q(a')<delta} 
    if(~C1)
        c=1;
        for i=1:m
            if((Q(i)-Q_primo)<delta)
                CS(c)=i;
                c=(c+1);
            end
        end
        disp('La condizione C1 non è verificata')
        disp('CS (compromise solution) = ');
        disp(CS);
    end
% SCATTER PLOT
    TP = figure();
    axis off
    hold on
    axis equal
    scatter(QS,QR, 40, 'green', 'fill');
% RIEMPIMENTO AREA DI INTERESSE PER LA SOLUZIONE OTTIMALE
    fill([2*Q_primo 0 0 (2*Q_primo+2*delta)],[0 2*Q_primo (2*Q_primo+2*delta) 0],'y');
% IDENTIFICAZIONE CITTA' IN CS
    for i=1:max(size(CS))
        scatter(QS(CS(i)),QR(CS(i)), 40, 'r', 'fill');
    end
% QUADRATO UNITARIO
    line( [0 0], [0 1], 'Color', 'black')
    line( [0 1], [1 1], 'Color', 'black')
    line( [1 1], [1 0], 'Color', 'black')
    line( [1 0], [0 0], 'Color', 'black')
% RETTA DI EQUAZIONE QS+QR=2*(Qmin)
   line([0 (2*Q_primo)], [(2*Q_primo) 0], 'Color', 'red') %retta pasante per a'
% RETTA DI EQUAZIONE QS+QR=2*(Qmin+delta)
   line([0 (2*Q_primo+2*delta)], [(2*Q_primo+2*delta) 0], 'Color', 'red')%retta interesse per C1
% RAPPRESENTAZIONE INDICI DELLE ALTERNATIVE PLOTTATE CON FUNZIONE scatter.
    for i = 1:m
        text( QS(i), QR(i), int2str(i), 'Color', 'black', 'FontSize', 12);
    end
% STAMPA DEL GRAFICO FINALE
    print( TP, 'TP',  '-dpng' );



