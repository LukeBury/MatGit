function qErr = quaternionError(qDes,qCur)
% Henri Kjellberg derived 
% Find the quaternion error (qErr) between desired quaternion (qDes) and current
% quaternion (qCur): q4 = scalar.

qDes = [2;3;4;3];
qCur = [5;7;3;4];


qErr = [qDes(4)  qDes(3) -qDes(2) -qDes(1);
       -qDes(3)  qDes(4)  qDes(1) -qDes(2);
        qDes(2) -qDes(1)  qDes(4) -qDes(3);
        qDes(1)  qDes(2)  qDes(3)  qDes(4)] * qCur;
    

qErr


end

    
    



