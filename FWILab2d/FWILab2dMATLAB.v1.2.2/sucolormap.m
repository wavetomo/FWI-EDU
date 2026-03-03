function sucolormap(ax,index)
% sucolormap(ax,index)
% set the axis 'ax' to a su hsv colormap (index>=0)
%   or conventional other colormap (index<0)
% index:
% -6 -> MapRed2White2Black 
% -5 -> MapRedYellowGreenBlueWhite (recommend for velocity display)
% -4 -> MapJetBlue2Green2Red
% -3 -> MapJetRed2Green2Blue (recommend for velocity display;
%                             this is from matlab's jet colormap)
% -2 -> MapBlack2White2Brown
% -1 -> MapBrown2White2Black (recommend for reflector display)
% 0 -> MapBlack2Gray2White
% 1 -> MapWhite2Gray2Black
% 2 -> MapBlue2Green2Red
% 3 -> MapRed2Green2Blue (recommend for velocity display)
% 4 -> MapBrown2Green2Blue
% 5 -> MapRed2White2Blue (recommend for residual display)
% 6 -> MapBlue2White2Red
% 7 -> MapWhite2Red2Yellow2Green2Blue
% 8 -> MapWhite2Green2Blue
% 9 -> MapBlue2Green2Yellow2Red2White
% 10 -> MapBlue2Green2White
% 11 -> MapBlue2White2Green
% 12 -> MapYellow2Red2Brown
% 13 -> MapBrown2Red2Yellow
% 14 -> MapRed2Yellow2Brown
if nargin<2
    index=0;
end
if nargin<1
    disp('the number of input parameters is one or two!');
    return;
end
if index>14
    index=mod(index,14);
end
if index<-6
    disp('index is out of range [-6 14]!');
    return;
end
if mod(index,1)~= 0
    disp('index must be an integer number!');
    return;
end
if index==-6
Map0=load('sucolormap.mat','MapRed2White2Black');
colormap(ax,Map0.MapRed2White2Black);
end
if index==-5
Map0=load('sucolormap.mat','MapRedYellowGreenBlueWhite');
colormap(ax,Map0.MapRedYellowGreenBlueWhite);
end
if index==-4
Map0=load('sucolormap.mat','MapJetBlue2Green2Red');
colormap(ax,Map0.MapJetBlue2Green2Red);
end
if index==-3
Map0=load('sucolormap.mat','MapJetRed2Green2Blue');
colormap(ax,Map0.MapJetRed2Green2Blue);
end
if index==-2
Map0=load('sucolormap.mat','MapBlack2White2Brown');
colormap(ax,Map0.MapBlack2White2Brown);
end
if index==-1
Map0=load('sucolormap.mat','MapBrown2White2Black');
colormap(ax,Map0.MapBrown2White2Black);
end
if index==0
Map0=load('sucolormap.mat','MapBlack2Gray2White');
colormap(ax,Map0.MapBlack2Gray2White);
end
if index==1
Map1=load('sucolormap.mat','MapWhite2Gray2Black');
colormap(ax,Map1.MapWhite2Gray2Black);
end
if index==2
Map2=load('sucolormap.mat','MapBlue2Green2Red');
colormap(ax,Map2.MapBlue2Green2Red);
end
if index==3
Map3=load('sucolormap.mat','MapRed2Green2Blue');
colormap(ax,Map3.MapRed2Green2Blue);
end
if index==4
Map4=load('sucolormap.mat','MapBrown2Green2Blue');
colormap(ax,Map4.MapBrown2Green2Blue);
end
if index==5
Map5=load('sucolormap.mat','MapRed2White2Blue');
colormap(ax,Map5.MapRed2White2Blue);
end
if index==6
Map6=load('sucolormap.mat','MapBlue2White2Red');
colormap(ax,Map6.MapBlue2White2Red);
end
if index==7
Map7=load('sucolormap.mat','MapWhite2Red2Yellow2Green2Blue');
colormap(ax,Map7.MapWhite2Red2Yellow2Green2Blue);
end
if index==8
Map8=load('sucolormap.mat','MapWhite2Green2Blue');
colormap(ax,Map8.MapWhite2Green2Blue);
end
if index==9
Map9=load('sucolormap.mat','MapBlue2Green2Yellow2Red2White');
colormap(ax,Map9.MapBlue2Green2Yellow2Red2White);
end
if index==10
Map10=load('sucolormap.mat','MapBlue2Green2White');
colormap(ax,Map10.MapBlue2Green2White);
end
if index==11
Map11=load('sucolormap.mat','MapBlue2White2Green');
colormap(ax,Map11.MapBlue2White2Green);
end
if index==12
Map12=load('sucolormap.mat','MapYellow2Red2Brown');
colormap(ax,Map12.MapYellow2Red2Brown);
end
if index==13
Map13=load('sucolormap.mat','MapBrown2Red2Yellow');
colormap(ax,Map13.MapBrown2Red2Yellow);
end
if index==14
Map14=load('sucolormap.mat','MapRed2Yellow2Brown');
colormap(ax,Map14.MapRed2Yellow2Brown);
end

    



