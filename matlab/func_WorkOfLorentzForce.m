function WLorentz = func_WorkOfLorentzForce(Bx,By,Bz,Jx,Jy,Jz,Vx,Vy,Vz)
    JxB_x = Jy.*Bz - Jz.*By;
    JxB_y = Jz.*Bx - Jx.*Bz;
    JxB_z = Jx.*By - Jy.*Bx;

    WLorentz = Vx.*JxB_x + Vy.*JxB_y + Vz.*JxB_z;
end