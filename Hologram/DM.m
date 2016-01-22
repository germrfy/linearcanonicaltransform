% Propagation over a short distance: direct method

%	d: distance of propagation
%	ui: wavefront before propagation
%   delx,dely: hologram sampling step

% Modified 2010-12-1

function[retVar]=DM(d,ui,delx,dely,lambda)
%  Defining input extent
[Nx,Ny]=size(ui);
Ax=Nx*delx;
Ay=Ny*dely;
Axvector=-Ax/2:delx:Ax/2-delx;
Ayvector=-Ay/2:dely:Ay/2-dely;
[AX,AY]=meshgrid(Axvector,Ayvector);

%  Defining spatial frequency extent
Bx=abs(lambda*d/delx);
By=abs(lambda*d/dely);
% Bx=lambda*d/delx;
% By=lambda*d/dely;
delBx=Bx/Nx;
delBy=By/Ny;
Bxvector=-Bx/2:delBx:Bx/2-delBx;
Byvector=-By/2:delBy:By/2-delBy;
[BX,BY]=meshgrid(Bxvector,Byvector);



     chirp=exp((1j*pi/(lambda*d)).*(AX.^2+AY.^2));
     const=exp(1j*pi/(lambda*d).*(BX.^2+BY.^2));
   
     retVar=const.*fftshift(fft2(fftshift(ui.*chirp)))/sqrt(Nx*Ny);


     if d<0
          retVar=rot90(retVar,2);
          retVar=circshift(retVar,[1,1]);
     end
     



 