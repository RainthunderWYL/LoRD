function out = Tool_ShowProjectedField(haxis,B1,B2,B3,x1,x2,x3,r0,r1,...
    MSP_ScaleRatio,NumRefine,MSPDispScale)

e3i = r1-r0;
e3i = e3i/norm(e3i);

e3.e1 = e3i(1);
e3.e2 = e3i(2);
e3.e3 = e3i(3);

[e1,e2,e3] = func_VectorLocalFrame(e3);

e1i = [e1.e1,e1.e2,e1.e3];
e2i = [e2.e1,e2.e2,e2.e3];

out = func_ProjectFieldProf(haxis,B1,B2,B3,x1,x2,x3,r0,e1i,e2i,e3i,...
    MSP_ScaleRatio,NumRefine,MSPDispScale);
end