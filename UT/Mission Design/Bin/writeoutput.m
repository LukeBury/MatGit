function writeoutput(dv,dvmag,dve,JD,tof,rsoi,vinfout,vinfin,vinfpl,rP,vP,...
    vPplhat,vd,va,vplhat,coet,dvRA,dvDEC,rm,rmmagr,vout,vin,voutmag,vinmag,rpout,rpin,vpout,...
    vpin,vpoutmag,vpinmag,rpfb,rpmag,rfbc,rpfbmagr,betaout,betain,delout,delin,deltaoutr,...
    deltainr,tsoiout,tsoiin,Tp,mupl,rpl,rap,dp,wp,mup,coep,coenp,rnb,rnb0,event,...
    optimization,npl,grade,rev,revig,vinfavmag,ilv,altp,optr,optnode,...
    dvopt,optoptions,ephtype,File,nodes,trajs,extime,vinfoutmag,vinfinmag,output,flag)

if optimization~=0
output=struct2cell(output);
end

file=strcat(File,'.xlsm');
copyfile('Bin/OutputTemp.xlsm',file)

Excel = actxserver ('Excel.Application'); 
Excel.Workbooks.Open(fullfile(pwd,strcat('\',file)));

% tic
xlswrite1(file,{file},1,'B1')
xlswrite1(file,{date},1,'D1')
xlswrite1(file,{extime},1,'F1')
xlswrite1(file,{flag},1,'H1')

xlswrite1(file,{cell2mat(output(1))},1,'B2')
xlswrite1(file,{cell2mat(output(2))},1,'D2')
xlswrite1(file,{cell2mat(output(3))},1,'F2')
xlswrite1(file,{cell2mat(output(4))},1,'H2')

if ephtype==1
    xlswrite1(file,{'DE405'},1,'B3')
elseif ephtype==2
    xlswrite1(file,{'PolyFit 3o'},1,'B3')
end
% toc

% tic 

    
%  xlswrite1(file,{strcat('Node',num2str(n))},1,strcat('A',num2str(2*n+3)))
xlswrite1(file,reshape(nodes,1,[]),'Raw Data',strcat('B',num2str(1)))
xlswrite1(file,reshape(dv,1,[]),'Raw Data',strcat('B',num2str(2)))
xlswrite1(file,reshape(dvmag,1,[]),'Raw Data',strcat('B',num2str(3)))
xlswrite1(file,reshape(dve,1,[]),'Raw Data',strcat('B',num2str(4)))
xlswrite1(file,reshape(JD,1,[]),'Raw Data',strcat('B',num2str(5)))
xlswrite1(file,reshape(tof,1,[]),'Raw Data',strcat('B',num2str(6)))
xlswrite1(file,reshape(rsoi,1,[]),'Raw Data',strcat('B',num2str(7)))
xlswrite1(file,reshape(vinfout,1,[]),'Raw Data',strcat('B',num2str(8)))
xlswrite1(file,reshape(vinfin,1,[]),'Raw Data',strcat('B',num2str(9)))
xlswrite1(file,reshape(vinfpl,1,[]),'Raw Data',strcat('B',num2str(10)))
xlswrite1(file,reshape(rP,1,[]),'Raw Data',strcat('B',num2str(11)))
xlswrite1(file,reshape(vP,1,[]),'Raw Data',strcat('B',num2str(12)))
xlswrite1(file,reshape(vPplhat,1,[]),'Raw Data',strcat('B',num2str(13)))
xlswrite1(file,reshape(vd,1,[]),'Raw Data',strcat('B',num2str(14)))
xlswrite1(file,reshape(va,1,[]),'Raw Data',strcat('B',num2str(15)))
xlswrite1(file,reshape(vplhat,1,[]),'Raw Data',strcat('B',num2str(16)))
xlswrite1(file,reshape(rm,1,[]),'Raw Data',strcat('B',num2str(17)))
xlswrite1(file,reshape(vout,1,[]),'Raw Data',strcat('B',num2str(18)))
xlswrite1(file,reshape(vin,1,[]),'Raw Data',strcat('B',num2str(19)))
xlswrite1(file,reshape(rpout,1,[]),'Raw Data',strcat('B',num2str(20)))
xlswrite1(file,reshape(rpin,1,[]),'Raw Data',strcat('B',num2str(21)))
xlswrite1(file,reshape(vpout,1,[]),'Raw Data',strcat('B',num2str(22)))
xlswrite1(file,reshape(vpin,1,[]),'Raw Data',strcat('B',num2str(23)))
xlswrite1(file,reshape(rpfb,1,[]),'Raw Data',strcat('B',num2str(24)))
xlswrite1(file,reshape(rpmag,1,[]),'Raw Data',strcat('B',num2str(25)))
xlswrite1(file,reshape(rfbc,1,[]),'Raw Data',strcat('B',num2str(26)))
xlswrite1(file,reshape(rpfbmagr,1,[]),'Raw Data',strcat('B',num2str(27)))
xlswrite1(file,reshape(betaout*180/pi,1,[]),'Raw Data',strcat('B',num2str(28)))
xlswrite1(file,reshape(betain*180/pi,1,[]),'Raw Data',strcat('B',num2str(29)))
xlswrite1(file,reshape(delout*180/pi,1,[]),'Raw Data',strcat('B',num2str(30)))
xlswrite1(file,reshape(delin*180/pi,1,[]),'Raw Data',strcat('B',num2str(31)))
xlswrite1(file,reshape(tsoiout,1,[]),'Raw Data',strcat('B',num2str(32)))
xlswrite1(file,reshape(tsoiin,1,[]),'Raw Data',strcat('B',num2str(33)))
xlswrite1(file,reshape(Tp,1,[]),'Raw Data',strcat('B',num2str(34)))
xlswrite1(file,reshape(mupl,1,[]),'Raw Data',strcat('B',num2str(35)))
xlswrite1(file,reshape(rpl,1,[]),'Raw Data',strcat('B',num2str(36)))
xlswrite1(file,reshape(rap,1,[]),'Raw Data',strcat('B',num2str(37)))
xlswrite1(file,reshape(dp,1,[]),'Raw Data',strcat('B',num2str(38)))
xlswrite1(file,reshape(wp,1,[]),'Raw Data',strcat('B',num2str(39)))
xlswrite1(file,reshape(mup,1,[]),'Raw Data',strcat('B',num2str(40)))
xlswrite1(file,reshape(coep,1,[]),'Raw Data',strcat('B',num2str(41)))
xlswrite1(file,reshape(coenp,1,[]),'Raw Data',strcat('B',num2str(42)))
xlswrite1(file,reshape(rnb,1,[]),'Raw Data',strcat('B',num2str(43)))
xlswrite1(file,reshape(event,1,[]),'Raw Data',strcat('B',num2str(44)))
xlswrite1(file,reshape(optimization,1,[]),'Raw Data',strcat('B',num2str(45)))
xlswrite1(file,reshape(npl,1,[]),'Raw Data',strcat('B',num2str(46)))
xlswrite1(file,reshape(grade,1,[]),'Raw Data',strcat('B',num2str(47)))
xlswrite1(file,reshape(rev,1,[]),'Raw Data',strcat('B',num2str(48)))
xlswrite1(file,reshape(revig,1,[]),'Raw Data',strcat('B',num2str(49)))
xlswrite1(file,reshape(vinfavmag,1,[]),'Raw Data',strcat('B',num2str(50)))
xlswrite1(file,reshape(ilv,1,[]),'Raw Data',strcat('B',num2str(51)))
xlswrite1(file,reshape(altp,1,[]),'Raw Data',strcat('B',num2str(52)))
xlswrite1(file,reshape(0,1,[]),'Raw Data',strcat('B',num2str(53)))
xlswrite1(file,reshape(optr,1,[]),'Raw Data',strcat('B',num2str(54)))
xlswrite1(file,reshape(optnode,1,[]),'Raw Data',strcat('B',num2str(55)))
xlswrite1(file,reshape(dvopt,1,[]),'Raw Data',strcat('B',num2str(56)))
xlswrite1(file,reshape(vinfoutmag,1,[]),'Raw Data',strcat('B',num2str(57)))
xlswrite1(file,reshape(vinfinmag,1,[]),'Raw Data',strcat('B',num2str(58)))
xlswrite1(file,reshape(rmmagr,1,[]),'Raw Data',strcat('B',num2str(59)))
xlswrite1(file,reshape(voutmag,1,[]),'Raw Data',strcat('B',num2str(60)))
xlswrite1(file,reshape(vinmag,1,[]),'Raw Data',strcat('B',num2str(61)))
xlswrite1(file,reshape(vpoutmag,1,[]),'Raw Data',strcat('B',num2str(62)))
xlswrite1(file,reshape(vpinmag,1,[]),'Raw Data',strcat('B',num2str(63)))
xlswrite1(file,reshape(deltaoutr,1,[]),'Raw Data',strcat('B',num2str(64)))
xlswrite1(file,reshape(deltainr,1,[]),'Raw Data',strcat('B',num2str(65)))
xlswrite1(file,reshape(coet,1,[]),'Raw Data',strcat('B',num2str(66)))
xlswrite1(file,reshape(dvRA,1,[]),'Raw Data',strcat('B',num2str(67)))
xlswrite1(file,reshape(dvDEC,1,[]),'Raw Data',strcat('B',num2str(68)))
xlswrite1(file,reshape(rnb0,1,[]),'Raw Data',strcat('B',num2str(69)))


Excel.Run('Config');
invoke(Excel.ActiveWorkbook,'Save'); 
Excel.Quit 
Excel.delete 
clear Excel

% toc

end