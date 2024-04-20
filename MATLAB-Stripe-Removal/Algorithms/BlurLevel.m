function blur = BlurLevel(original)

I = double(original);
[y x] = size(I);

% image filtering both sides
Hvh= ones(9,9)/81;
B_Ver_Hor = imfilter(I,Hvh);

 % Input image gradient new
 [Gmag, Gdir] = imgradient(I,'intermediate');
  Gmag = abs(Gmag);
  
 % Blur image gradient new
  [GmagB, GdirB] = imgradient(B_Ver_Hor,'intermediate');
  GmagB = abs(GmagB);
  
 % New Difference between blur and input
 T_Ver_Hor = Gmag - GmagB;
 
 % New   D from the paper
  V_Ver_Hor =max(0,T_Ver_Hor);
  
 % Sum pixels gradient original  % sf
  S_D_Ver_Hor = sum(sum(Gmag(2:y-1,2:x-1)));
  
 % Sum pixels gradient differences oirginal blurred %sd
 S_V_Ver_Hor = sum(sum(V_Ver_Hor(2:y-1,2:x-1)));
 
 % Result
 blur = (S_D_Ver_Hor-S_V_Ver_Hor)/S_D_Ver_Hor;
 