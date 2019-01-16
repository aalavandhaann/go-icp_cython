from py_goicp import GoICP, POINT3D, ROTNODE, TRANSNODE;
import numpy as np;
import time;

def loadPointCloud(filename):
    pcloud = np.loadtxt(filename, skiprows=1);
    plist = pcloud.tolist();
    p3dlist = [];
    for x,y,z in plist:
        pt = POINT3D(x,y,z);
        p3dlist.append(pt);
    return pcloud.shape[0], p3dlist;

goicp = GoICP();
rNode = ROTNODE();
tNode = TRANSNODE();
 
rNode.a = -3.1416;
rNode.b = -3.1416;
rNode.c = -3.1416;
rNode.w = 6.2832;
 
tNode.x = -0.5;
tNode.y = -0.5;
tNode.z = -0.5;
tNode.w = 1.0;

goicp.MSEThresh = 0.001;
goicp.trimFraction = 0.0;
  
if(goicp.trimFraction < 0.001):
    goicp.doTrim = False;

a_points = [POINT3D(0.0, 0.0, 0.0), POINT3D(0.5, 1.0, 0.0), POINT3D(1.0, 0.0, 0.0)];
b_points = [POINT3D(0.0, 0.0, 0.0), POINT3D(1.0, 0.5, 0.0), POINT3D(1.0, -0.5, 0.0)];

Nm, a_points = loadPointCloud('./tests/model_bunny.txt');
Nd, b_points = loadPointCloud('./tests/data_bunny.txt');

goicp.loadModelAndData(Nm, a_points, Nd, b_points);
#LESS DT Size = LESS TIME CONSUMPTION = HIGHER ERROR
goicp.setDTSizeAndFactor(10, 2.0);
goicp.setInitNodeRot(rNode);
goicp.setInitNodeTrans(tNode);

start = time.time();
print("Building Distance Transform...");
goicp.BuildDT();
print("REGISTERING....");
goicp.Register();
end = time.time();
print('TOTAL TIME : ', (end-start));
optR = np.array(goicp.optimalRotation());
optT = goicp.optimalTranslation();
optT.append(1.0);
optT = np.array(optT);

transform = np.empty((4,4));
transform[:3, :3] = optR;
transform[:,3] = optT;

print(optR);
print(optT);
print(transform);




