%--------------------------------------------------------------
% TOPSIS Method 
%--------------------------------------------------------------

clear
%------------------- Read in data
X = dlmread('Mglobal.txt');
disp('Dimensions: (m alternatives, n criteria)')
[m,n] = size(X)
disp('Criteria directions (1: max, -1: min)')
Dir = dlmread('Dir.txt')
disp('Criteria weights ')
W = dlmread('W.txt')
w = W/norm(W,1);

%--------------------------------------------------------------
% Normalization: mapping onto [0,1]
% NOTE: not used, overwritten below!
%--------------------------------------------------------------
for j = 1:n
  val_min = min(X(:,j));
  val_max = max(X(:,j));
  for i = 1:m
    R(i,j) = (X(i,j)-val_min)/(val_max-val_min);
  end 
end

%--------------------------------------------------------------
% Normalization: division by 2-norm
for j = 1:n
   R(:,j) = X(:,j) / norm( X(:,j), 2 );
end

%--------------------------------------------------------------
% Weighted matrix V
for j=1:n
  V(:,j) = w(j) * R(:,j);
end

%--------------------------------------------------------------
% Compute Ideal and ANTI-Ideal
for j=1:n
  if Dir(j)==1
     I(j)=max(V(:,j));
     AI(j)=min(V(:,j));
  else
     I(j)=min(V(:,j));
     AI(j)=max(V(:,j));
  end
end

%--------------------------------------------------------------
% Compute Euclidean distances from Ideal and Anti-Ideal
% Compute SI (distance from ideal I) SAI (distance from anti-ideal AI)
for i=1:m
  SI(i) = norm(I-V(i,:),2); 
  SAI(i) = norm(V(i,:)-AI,2);
end

%--------------------------------------------------------------
% Compute Index C
for i=1:m
  C(i)= SAI(i)/(SAI(i)+SI(i));
end

%--------------------------------------------------------------
% Sort index values and print ranking
disp('TOPSIS Ranking (decreasing index C)')
Ranking = sortrows( [ C' (1:m)' ] , 'descend' );
for i = 1:m
    st = sprintf("%2d) \t%2d \t%f", i, Ranking(i,2), Ranking(i,1) );
    disp(st)
end


%--------------------------------------------------------------
% cmpute normalized distances from Ideal and Anti-Ideal 
L = norm(I-AI,2);
NSI = SI/L;
NSAI = SAI/L;

%--------------------------------------------------------------
% Plot alternatives on the (NSI,NSAI) plane
TP = figure();
hold on
axis equal
axis off
scatter(NSI,NSAI, 10, 'b', 'fill');
% draw diagonal
plot([0 1], [1 0], 'Color', 'r')
% draw bounding box
% Draw bounding box and axes
line( [0 0], [0 1], 'Color', 'black')
line( [0 1], [1 1], 'Color', 'black')
line( [1 1], [1 0], 'Color', 'black')
line( [1 0], [0 0], 'Color', 'black')
% draw TOPSIS iso-index ines, with corresponding alternative index
for i = 1:m
    co = [NSI(i),NSAI(i)];
    co = sqrt(2) * co/norm(co,2);
    line([0 co(1)], [0, co(2)], 'color', 'green' );
    text( co(1), co(2), int2str(i), 'Color', 'b', 'FontSize', 8  );
end
% plot circle containing points
angle = 0:pi/100:pi/2;
x = cos(angle);
y = sin(angle);
plot(x,y, 'Color', 'r')
% save picture
print( TP, 'TP',  '-dpng' );

%--------------------------------------------------------------





