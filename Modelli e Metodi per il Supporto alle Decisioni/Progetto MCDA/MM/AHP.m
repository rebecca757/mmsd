% *** Metodologia AHP       (Versione 15/10/23) ***
% *** INPUT matrice M quadrata, di interi:      ***
% ***       i reciproci sono sostituiti da      ***
% ***       0 e inseriti dal programma          ***
% *** OUTPUT autovettore princiale (modulo 1);  ***
% ***        vettore w di pesi, CI; RI se n<10  ***
%
% Valori CRI (dalla letteratura)
CRI(9) = 1.45;
CRI(8) = 1.41;
CRI(7) = 1.32;
CRI(6) = 1.24;
CRI(5) = 1.12;
CRI(4) = 0.9;
CRI(3) = 0.58;
CRI(2) = 0;
CRI(1) = 0;
%
% Lettura da file di una matrice quadrata (scritta per righe)
% Si assume che la matrice contenga solo valori interi: i valori
% frazionari sono sostituiti da zero, e sono ricalcolati, cioè,
% se M(i,j) == 0 si pone M(i,j) = 1 / M(j,i) 
% I valori sulla diagonale vengono comunque fissati a 1
disp(' ')
disp('Lettura matrice M')
M = dlmread('A.txt');
[m,n] = size(M);
disp('Dimensione:')
n
%
% Si trasforma M in una matrice reciproca
for i = 1:n 
    M(i,i) = 1;
    for j = (1):(i-1)
        if M(i,j) == 0
           if M(j,i) == 0
               disp('ERRORE: matrice non corretta')
               return
           end
           M(i,j) = 1 / M(j,i);
        else
           M(j,i) = 1 / M(i,j);
        end
    end
end
disp('Matrice reciproca (approssimata)')
M
% Calcolo autovalore principale (colonna D) e rispettivo
% autovettore (D) 
disp('(Calcolo autovalore e autovettore)')
SV = ones(n,1);
[V,D] = eigs(M, 1, 'lm', 'StartVector', SV );
% chiamo l'autovalore principale
W = V;
if ( ones(1,n)*W < 0) 
    disp(' ')
    disp('Autovettore negativo: verso cambiato!')
    W = -1*W;
end
disp(' ')
disp('Autovettore massimo W (modulo 1)')
disp(W')
% normalizzazione: somma dei valori = 1
disp('Vettore dei pesi w (somma 1)')
w = W / norm(W,1);
disp(w')
%
% Calcolo CI e RI
disp('')
disp('Autovalore principale (lambda1)')
lambda1 = D(1,1)
disp('')
disp('Indice di Inconsistenza:')
CI = ( lambda1 - n ) / (n-1)
if n < 10
    disp('')
    disp('Rapporto di Inconsistenza:')
    RI = CI / CRI(n)
else 
    disp('')
    disp('Rapporto di Inconsistenza non calcolabile')
end
disp('')
   