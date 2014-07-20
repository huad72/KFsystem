#include "testApp.h"
#define KINECT_WIDTH 640
#define KINECT_HEIGHT 480

#define SHORT_DMAX 32767

#define SIZE_TRIDE 0.015625//0.015625
#define TSDF_SIZE 256 //128
#define WRONG_DATA 3

#define THRESHOLD_D 0.1 // 要修改
#define THRESHOLD_A 20 // 角度，要修改
//ofMatrix3x3 getRotationMatrix(const ofMatrix4x4& Tmaxtrix);
//--------------------------------------------------------------
void testApp::testh()
{
	//ofMatrix4x4 vec=ofMatrix4x4(2,0,0,0,0,3,0,0,0,0,4,0,0,0,0,5);
	//cout<<vec<<endl;
	//vec=vec.getInverse();
	//ofVec4f is=ofVec4f(1.0,1.0,1.0,1.0);
	//ofVec4f it = vec*is;
	//cout<<vec<<endl;

	//union name{
	//	uchar c[4];
	//	float a;
	//};
	//name na;
	//na.a=1.2;
 ////   na.c[0]=255;
	////na.c[1]=111;
	////na.c[2]=159;
	////na.c[3]=63;
	//
	//cout<<na.a<<endl;
	//for(int i=0;i<4;i++)
	//{
	//	cout<<int(na.c[i])<<endl;
	//}

	//int a,b;
	//ofstream out;
	//out.open("pointdata1.txt",ios::trunc);
	//for(int i=0;i<10;i++)
	//{
	//	//float k = 1.999*i;
	//	ofVec3f poi = ofVec3f(i,i*1.999,10);
	//	out<<poi.x<<" "<<poi.y<<" "<<poi.z<<" ";
	//}
	//out.close();
	//ifstream fin;
	//fin.open("pointdata1.txt");
	//for(int i=0;i<10;i++)
	//{
	//	float k[3];
	//	for(int j=0;j<3;j++)
	//	{
	//	fin>>k[j];
	//	}
	//	ofVec3f poi(k[0],k[1],k[2]);
	//	//tsdfData[i] = k;
	//	cout<<poi<<endl;
	//	//cout<<k<<" ";
	//}

	//freopen("data.txt","r",stdin);
	//for(int i=0;i<TSDF_SIZE*TSDF_SIZE*TSDF_SIZE;i++)
	//{
	//	//float k;
	//	scanf("%f",&tsdfData[i]);
	//	//tsdfData[i] = k;
	//	cout<<tsdfData[i]<<" ";
	//}
	//fclose(stdin);

	//cout<< numeric_limits<float>::quiet_NaN ()<<endl;
	//ofMatrix4x4 tk=ofMatrix4x4(1,2,3,4,5,-1,6,7,8,9,10,11,12,13,14,15);
	//cout<<tk<<endl;
	//ofMatrix3x3 ss = getRotationMatrix(tk);
	//cout<<ss<<endl;

	ofVec3f normal1 = ofVec3f(1,2,3);
	ofVec3f normal2 = ofVec3f(1,0,0);
	//float angle = normal1.angle(normal2);
	//cout<<angle<<endl;
	//compute_icp(5,normal1,normal2);
	//Eigen::MatrixXf startm = Eigen::MatrixXf::Zero(6,6);
	//cout<<startm<<endl;
	//ofMatrix4x4 cha = ofMatrix4x4(1,0,0,1,0,-1,0,2,0,0,1,3,0,0,0,1);
	//normal2 = cha * normal1;
	//cout<<normal2<<endl;
	//Eigen::Matrix<float,3,6> Gx;
	//Gx<<(Eigen::Matrix3f()<<0,-1,2,1,0,-3,-2,3,0).finished(),
	//	 Eigen::MatrixXf::Identity(3,3);
	//cout<<Gx<<endl;
}
//--------------------------文件读取----------------------------
void testApp::readglobalfile(char *filename,float *data)
{
	//freopen("data.txt","r",stdin);
	//for(int i=0;i<TSDF_SIZE*TSDF_SIZE*TSDF_SIZE;i++)
	//{
	//	//float k;
	//	scanf("%f",&tsdfData[i]);
	//	//tsdfData[i] = k;
	//	cout<<tsdfData[i]<<" ";
	//}
	//fclose(stdin);

	// 快速读取
	const int maxs = 6*1024*1024;
	char *buf = new char[maxs];
	freopen(filename,"rb",stdin);
	int len = fread(buf,1,maxs,stdin);
	buf[len] = '\0';
	float k = 0;
	bool neg = false;
	bool vig = false;
	float dec = 0.1;
	int i =0;
	int m = 1;
	for(char* p=buf;*p&&p-buf<len;p++)
	{
		if(*p == ' ')
		{
			if(neg)
				k*=-1;
			data[i]=k;
			i++;
			//cout<<k<<" ";
			neg = false;
			vig = false;
			dec = 0.1;
			k = 0;
		}
		else
		{
			if(vig)
			{
				k+=(*p-'0')*dec;
				dec*=0.1;
			}
			else
			{
				if(*p=='-')
					neg=true;
				if(*p=='.')
					vig=true;
				else
				{
					//k *= m;
					k =k*10 + (*p-'0');
				}
			}

		}
	}
	//cout<<i<<endl;
}
void readFileToArry(ofVec3f *pointsmap,char* filename)
{
	int realnum = 0;// 统计数量
	const int maxs = 6*1024*1024;
	char *buf = new char[maxs];
	freopen(filename,"rb",stdin);
	int len = fread(buf,1,maxs,stdin);
	buf[len] = '\0';
	float k;
	bool neg = false;
	bool vig = false;
	float dec = 0.1;
	int i = 0;
	int j = 0;
	float num_vector[3];
	for(char* p=buf;*p&&p-buf<len;p++)
	{
		if(*p == ' ')
		{
			if(neg)
				k*=-1;
			num_vector[j%3] = k;
			if((j+1)%3==0)
			{
				pointsmap[i].x=num_vector[0];
				pointsmap[i].y=num_vector[1];
				pointsmap[i].z=num_vector[2];
				//pointsmap[i].x=num_vector[0]*100;
				//pointsmap[i].y=num_vector[1]*100;
				//pointsmap[i].z=num_vector[2]*100;
				if(pointsmap[i].z>0)
					++ realnum;
				//cout<<pointsmap[i]<<endl;
				i++;
			}
			++j;
			//cout<<k<<" ";
			neg = false;
			vig = false;
			dec = 0.1;
		}
		else
		{
			if(vig)
			{
				k+=(*p-'0')*dec;
				dec*=0.1;
			}
			else
			{
			if(*p=='-')
				neg=true;
			if(*p=='.')
				vig=true;
			else
				k=*p-'0';
			}

		}
	}
	cout<<"no zero point num = "<<realnum<<endl;
}
// icp测试
void testTheIcp(const ofVec3f* pointsmap,ofVec3f* pointsmap_f)
{
	ofMatrix4x4 tm = ofMatrix4x4(0.98481,-0.17365,0,0,0.17365,0.98481,0,0,0,0,1,0,0,0,0,1);
	//ofMatrix3x3 tm = ofMatrix3x3(0.98481,-0.17365,0,0.17365,0.98481,0,0,0,1);
	for(int i = 0;i < KINECT_HEIGHT;  ++ i)
		for(int j = 0;j < KINECT_WIDTH; ++ j)
		{
			int num = i * KINECT_WIDTH + j;
			pointsmap_f[num] = pointsmap[num];
			if(pointsmap[num].z > 0)
			{
				ofVec4f newp = ofVec4f(pointsmap[num].x,pointsmap[num].y,pointsmap[num].z,1);
				ofVec4f newn = tm*newp;
				pointsmap_f[num].x = newn.x;
				pointsmap_f[num].y = newn.y;
				pointsmap_f[num].z = newn.z;
				//pointsmap_f[num].y += 0.1;
				//cout<<i<<endl;
			}
		}
}
//--------------------------------------------------------------
void testApp::setup() {
	testh();
	ofSetLogLevel(OF_LOG_VERBOSE);
	int count = 0;
	m_nui = NULL;
	bKinectInitSucessful = false;
	HRESULT hr=NuiGetSensorCount(&count);
	if(count<=0)
	{
		std::cout<<"No kinect sensor was found!!"<<endl;
		goto FINAL;
	}
	hr = NuiCreateSensorByIndex(0,&m_nui);
	if(FAILED(hr))
	{
		std::cout<<"create failed"<<endl; 
		goto FINAL;
	}
	hr=m_nui->NuiInitialize(NUI_INITIALIZE_FLAG_USES_DEPTH);
	//hr=m_nui->NuiInitialize(NUI_INITIALIZE_FLAG_USES_DEPTH_AND_PLAYER_INDEX);//NUI_INITIALIZE_FLAG_USES_DEPTH);
	if(FAILED(hr))
	{
		std::cout<<"NuiInitialize failed"<<endl; 
		goto FINAL;
	}
	
	depthStreamHandle=NULL;
	hr=m_nui->NuiImageStreamOpen(NUI_IMAGE_TYPE_DEPTH,NUI_IMAGE_RESOLUTION_640x480,0,2,NULL,&depthStreamHandle);
	//hr=m_nui->NuiImageStreamOpen(NUI_IMAGE_TYPE_DEPTH_AND_PLAYER_INDEX,NUI_IMAGE_RESOLUTION_320x240,0,2,NULL,&depthStreamHandle);
	if(FAILED(hr))
	{
		std::cout<<"Could not open color image stream video"<<endl; 
		goto FINAL;
	}

	//colorImg.allocate(KINECT_WIDTH, KINECT_HEIGHT);
	//grayImage.allocate(KINECT_WIDTH, KINECT_HEIGHT);
	//grayThreshNear.allocate(KINECT_WIDTH, KINECT_HEIGHT);
	//grayThreshFar.allocate(KINECT_WIDTH, KINECT_HEIGHT);
	
	nearThreshold = 255;
	farThreshold = 220;
	bThreshWithOpenCV = true;
	
	ofSetFrameRate(60);
	
	// 加入cvwindow
	cv::namedWindow("orginal");
	//cv::namedWindow("filter");
	//cv::namedWindow("downonce");
	//cv::namedWindow("downtwice");

	// start from the front
	bDrawPointCloud = false;
	tsdfFirstTime = true;

	// 我的改变
	test=true;
	test2=false;// 实验区控制
	test3=true;
	draww=false;
	countnum = 0;
	photonum = 0;
	// 光线投影点云数量
	finalPointNum = 0;

	// 摄像机内参
	the_K_cam=ofMatrix3x3(589.3,0,323,0,590.4,244.6,0,0,1);
	invKcam = the_K_cam.inverse(the_K_cam);
	//the_K_cam.invert();
	std::cout<<invKcam<<endl;

	//std::cout<<the_K_cam<<std::endl;
	//m_points=(ofVec3f*)malloc(sizeof(ofVec3f)*245*325);
	
	//ofVec3f vec[3][3];
	//int n=3;
	//for(int i=0;i<n;i++)
	//{
	//	for(int j=0;j<n;j++)
	//	{
	//		vec[i][j]=ofVec3f(i,j,i);
	//	}
	//}
	//compute_normal(&vec[0][0],3,3);

	depth_float= new float[KINECT_WIDTH*KINECT_HEIGHT];

	// 读取存在文件中的tsdf数据存入内存（注：方便调试用，以后舍弃）
	//tsdfData = new float[TSDF_SIZE*TSDF_SIZE*TSDF_SIZE];
	tsdfWeightData = new short2[TSDF_SIZE*TSDF_SIZE*TSDF_SIZE];

	// ICP得到的变换矩阵
	final_A = Eigen::MatrixXf::Zero(6,6);
	final_B = Eigen::MatrixXf::Zero(6,1);

	clock_t t1=clock();
	//readglobalfile("data.txt",tsdfData);
	//readglobalfile("0.txt",depth_float);
	//readFileToArry(pointsmap_orignal,"testdata1.txt");
	//readFileToArry(pointsmap_final,"testdata1.txt");
	//testTheIcp(pointsmap_orignal,pointsmap_final);
	clock_t t2=clock();
	double ti=(double)(t2-t1)/CLOCKS_PER_SEC;
	//cout<<ti<<endl;

	bKinectInitSucessful = true;

	// 摄像机位置控制
	x_back = 0;
	y_back = 0;
	z_back = -650;
FINAL:
	if(FAILED(hr))
	{
		std::cout<<"待做撤销工作"<<endl;
		if(m_nui != NULL)
		{
			m_nui->NuiShutdown();
			m_nui->Release();
			m_nui = NULL;
		}
	}
}
//----------------------矩阵乘法-废弃---------------------------
//ofVec3f testApp::nmmul(const ofMatrix3x3 &mat1,const ofVec3f &mat2) const{
//
//	ofVec3f resultmat(0);
//	resultmat.x=mat1.a*mat2.x+mat1.b*mat2.y+mat1.c*mat2.z;
//	resultmat.y=mat1.d*mat2.x+mat1.e*mat2.y+mat1.f*mat2.z;
//	resultmat.z=mat1.g*mat2.x+mat1.h*mat2.y+mat1.i*mat2.z;
//	return resultmat;
//}
//-----------------------出界判断--------------------------------
bool isvalid(const ofVec3f &gPostion)
{
	if(gPostion.x<0||gPostion.x>TSDF_SIZE-1)
	{
		//cout<<"x ="<<gPostion.x<<" is wrong"<<endl;
		return false;
	}
	if(gPostion.y<0||gPostion.y>TSDF_SIZE-1)
	{
		//cout<<"y ="<<gPostion.y<<" is wrong"<<endl;
		return false;
	}
	if(gPostion.z<0||gPostion.z>TSDF_SIZE-1)
	{
		//cout<<"z ="<<gPostion.z<<" is wrong"<<endl;
		return false;
	}
	return true;
}
//-----------------------tsdf数据存取----------------------------
void testApp::pack_tsdf(const float &tsdf,const int &weight,const int num)
{
	if(num<0&&num>=TSDF_SIZE*TSDF_SIZE*TSDF_SIZE)
	{
		cout<<"pack tsdf:num is wrong"<<endl;
	}
	else
	{
		int mixtsdf = max(-SHORT_DMAX,min(SHORT_DMAX,(int)(tsdf*SHORT_DMAX)));
		tsdfWeightData[num].tsdf = mixtsdf;
		tsdfWeightData[num].weight = weight;
	}
}
void testApp::unpack_tsdf(float &tsdf,int &weight,const int num)
{
	if(num<0&&num>=TSDF_SIZE*TSDF_SIZE*TSDF_SIZE)
	{
		cout<<"pack tsdf:num is wrong"<<endl;
	}
	else
	{
		tsdf = (float)tsdfWeightData[num].tsdf / SHORT_DMAX;
		weight = tsdfWeightData[num].weight;
	}
}
//-----------------------体元位置变换----------------------------
ofVec3f testApp::getTheGlobalPostion(const ofVec3f &gPoint)
{
	int x,y,z;
	x = (int)(gPoint.x/SIZE_TRIDE+0.5);
	y = (int)(gPoint.y/SIZE_TRIDE+0.5);
	z = (int)(gPoint.z/SIZE_TRIDE+0.5);

	//x =min(max((int)0,x),TSDF_SIZE-1);
	//y =min(max((int)0,y),TSDF_SIZE-1);
	//z =min(max((int)0,z),TSDF_SIZE-1);
	//cout<<x<<" "<<y<<" "<<z<<endl;
	return ofVec3f(x,y,z);
}
//-------------------------读取体元数据--------------------------
float testApp::getTsdfData(const int x,const int y,const int z)
{
	if(x<0||x>TSDF_SIZE-1)
	{
		std::cout<<"x<0"<<endl;
		return WRONG_DATA;
	}
	if(y<0||y>TSDF_SIZE-1)
	{
		std::cout<<"y<0"<<endl;
		return WRONG_DATA;
	}
	if(z<0||z>TSDF_SIZE-1)
	{
		std::cout<<"z<0"<<endl;
		return WRONG_DATA;
	}
	int num = x+y*TSDF_SIZE+z*TSDF_SIZE*TSDF_SIZE;
    float tsdf; //= tsdfData[num];
	int weight;
	unpack_tsdf(tsdf,weight,num);
	return tsdf;
}
//-------------------------三线性插值----------------------------
float testApp::triLinary(const ofVec3f &point)
{
	ofVec3f g = getTheGlobalPostion(point);

	if(g.x<=0||g.x>=TSDF_SIZE-1)
		return numeric_limits<float>::quiet_NaN();
	if(g.y<=0||g.y>=TSDF_SIZE-1)
		return numeric_limits<float>::quiet_NaN();
	if(g.z<=0||g.z>=TSDF_SIZE-1)
		return numeric_limits<float>::quiet_NaN();

	float vx = (g.x)*SIZE_TRIDE;
	float vy = (g.y)*SIZE_TRIDE;
	float vz = (g.z)*SIZE_TRIDE;
	g.x = (vx > point.x)? g.x-1:g.x;
	g.y = (vy > point.y)? g.y-1:g.y;
	g.z = (vz > point.z)? g.z-1:g.z;

	float a = (point.x-(g.x*SIZE_TRIDE))/SIZE_TRIDE;
	float b = (point.y-(g.y*SIZE_TRIDE))/SIZE_TRIDE;
	float c = (point.z-(g.z*SIZE_TRIDE))/SIZE_TRIDE;
	//cout<<a<<"　"<<b<<" "<<c<<endl;
	float f = getTsdfData(g.x,g.y,g.z) * (1 - a) * (1 - b) * (1 - c) +
		      getTsdfData(g.x + 1,g.y,g.z) * a * (1 - b) * (1 - c) +
		      getTsdfData(g.x,g.y + 1,g.z) * (1 - a) * b * (1 - c) +
		      getTsdfData(g.x,g.y,g.z + 1) * (1 - a) * (1 - b) * c +
		      getTsdfData(g.x + 1,g.y + 1,g.z) * a * b * (1 - c) +
		      getTsdfData(g.x + 1,g.y,g.z + 1) * a * (1 - b) * c +
		      getTsdfData(g.x,g.y + 1,g.z + 1) * (1 - a) * b * c +
		      getTsdfData(g.x + 1,g.y + 1,g.z + 1) * a * b * c;
	//cout<<getTsdfData(g.x,g.y,g.z)<<" "<<getTsdfData(g.x+1,g.y,g.z)<<" "<<getTsdfData(g.x,g.y+1,g.z)<<endl;
	//cout<<getTsdfData(g.x,g.y,g.z+1)<<"　"<<getTsdfData(g.x+1,g.y+1,g.z)<<" "<<getTsdfData(g.x+1,g.y,g.z+1)<<endl;
	//cout<<getTsdfData(g.x,g.y+1,g.z+1)<<"　"<<getTsdfData(g.x+1,g.y+1,g.z+1)<<endl;
	return f;
}
//-----------------------类型转换-废弃---------------------------
void testApp::changeType(const cv::Mat &depthimage,const int rows,const int cols)
{
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<cols;j++)
		{
			//depth_float[i*cols+j]=(float)depthimage.ptr<uchar>(i)[j];
			if(depth_float[i*cols+j]>0)
			{
				cout<<depth_float[i*cols+j]<<" ";
			}
		}
	}
}
//-------------------------法向计算------------------------------
void testApp::compute_normal(const ofVec3f *vec,const int rows,const int cols,ofVec3f *normalmap)
{
	ofVec3f normal1(0),normal2(0);
	const ofVec3f *line_one=vec;
	const ofVec3f *line_two=vec+cols;
	
	clock_t t1=clock();
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<cols;j++)
		{
			if(j!=cols-1)
			{
				normal1=line_one[j+1]-line_one[j];
				normal2=line_two[j]-line_one[j];
			}
			else if(j==cols-1)
			{
				normal1=line_one[j]-line_one[j-1];
				normal2=line_two[j]-line_one[j];
			}
			//cout<<normal1<<" "<<normal2<<endl;
			int num=i*cols+j;
			if(vec[num].z>0)
			{
				normalmap[num]=normal2.cross(normal1);
				normalmap[num]=normalmap[num].normalized();
			}
			else
			{
				normalmap[num]=ofVec3f(0.0,0.0,0.0);
			}
			//if(normalmap[num].z!=0)
			//{
			//	cout<<normalmap[num]<<endl;
			//}
		}
		if(i!=rows-2)
		{
			line_two+=cols;
			line_one+=cols;
		}
	}
	clock_t t2=clock();
	double ti=(double)(t2-t1)/CLOCKS_PER_SEC;
	//cout<<ti<<endl;
}
//--------------------------顶点计算-----------------------------
void testApp::compute_points(const cv::Mat &depthimage,const int rows,const int cols,ofVec3f* pointsmap)
{
	// 存入文件
	//ofstream out;
	//out.open("testdata6.txt",ios::trunc);

	ofVec3f res(0);
	ofVec3f pos(0,0,1);
	float num=0;
	int num2=0;
	float maxnum=0;;
	float minnum=4096.0;

	int nums=0;

	clock_t t1=clock();//开始计时
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<cols;j++)
		{
			pos.x=j;
			pos.y=i;
			//num2=i*cols+j;
			num=depth_float[num2];
			if(num>0)
			{
				nums++;
				pointsmap[num2]=nmmul(invKcam,pos);
				pointsmap[num2]*=num;
				//if(i==300)
				//{
				//	cout<<i<< " "<<j<<" "<<pointsmap[num2]<<endl;
				//}
			}
			else
			{
				pointsmap[num2]=ofVec3f(0.0,0.0,0.0);
			}
			// 数据写入文件
			//if(pointsmap[num2].z<1.8)
			//	out<<pointsmap[num2].x<<" "<<pointsmap[num2].y<<" "<<pointsmap[num2].z<<" ";
			//else
			//{
			//	out<<0<<" "<<0<<" "<<0<<" ";
			//}

			num2++;
			if(num>0)
			{
				if(num>maxnum)
					maxnum=num;
				if(num<minnum)
					minnum=num;
				//cout<<i<<" "<<j<<" "<<num<<endl;
			}
			if(rows==KINECT_HEIGHT)
			{
				if(num>0)
					isvalid_orignal[num2]=true;
				else
					isvalid_orignal[num2]=false;
			}
		}
	}
	clock_t t2=clock();
	double ti=(double)(t2-t1)/CLOCKS_PER_SEC;
	//cout<<nums<<endl;
	//cout<<ti<<endl;
	//cout<<maxnum<<" "<<minnum<<" "<<maxnum-minnum<<endl;

	//out.close();

}
//--------------------------tsdf计算-----------------------------
void testApp::compute_tsdf(const ofMatrix4x4 &t_g,cv::Mat &depthimage,const int g_size,const ofVec4f &camP)
{
	//ofVec3f iy = nmmul(invKcam,ofVec3f(0,0,1));
	
	// 写入文件
	//ofstream out;
	//out.open("data.txt",ios::trunc);

	cv::Mat color(g_size,g_size,CV_8UC1);
	cv::Mat bigcol;//(KINECT_HEIGHT,KINECT_WIDTH,CV_8UC1);
	cv::Mat bbcol;
	ofMatrix4x4 invTrans = t_g.getInverse();
	float changenum = 0.0;
	float dist = 0.0;
	float tsdf = 0.0;
	int wight = 0;
	float pretsdf = 0.0;
	int preweight = 0;
	int i,j,k;
	int num = 0;
	for(k=0;k<g_size;k++)
	{
	for(j=0;j<g_size;j++)
	for(i=0;i<g_size;i++)
	{
		//ofVec4f spacePostion = ofVec4f(2.0,2.0,1,1);
		ofVec4f spacePostion = ofVec4f(-2.0+i*SIZE_TRIDE,-2.0+j*SIZE_TRIDE,0.5+k*SIZE_TRIDE,1.0);
		
		ofVec4f camPostion = invTrans*spacePostion;
		ofVec3f dcamPostion = ofVec3f(camPostion.x,camPostion.y,camPostion.z);
		ofVec3f picPosition = nmmul(the_K_cam,dcamPostion);
		ofVec2f depthPostion = ofVec2f(picPosition.x/picPosition.z,picPosition.y/picPosition.z);
		float minus = 0;

		int x = (int)(max(minus,depthPostion.x)+0.5);
		int y = (int)(max(minus,depthPostion.y)+0.5);
		//cout<<x<<" "<<y<<endl;

		x = min(x,KINECT_WIDTH-1);
		y = min(y,KINECT_HEIGHT-1);
		//cout<<x<<" "<<y<<endl;

		ofVec3f rdepthPostion = ofVec3f(x,y,1);
		changenum = nmmul(invKcam,rdepthPostion).length();
		dist = spacePostion.distance(camP);
		float dd= depth_float[x+y*KINECT_WIDTH];
		if(dd>0)
			tsdf =-(dist/changenum-depth_float[x+y*KINECT_WIDTH]);
		else
			tsdf = 100;
		if(tsdf>0)
		{
			tsdf = min(1.0,tsdf/(SIZE_TRIDE*2));
		}
		else
		{
			tsdf = max(-1.0,tsdf/(SIZE_TRIDE*2));
		}
		// 写入文件与存入内存
		//out<<tsdf<<" ";
		//tsdfData[num] = tsdf;
		num++;
		if(k==g_size/2+100)
		{
			uchar *ptr = color.ptr<uchar>(j);
			uchar col = 0;
			if(tsdf > 0)
			{
				if(tsdf==1.0)
					col = 255;
				else
					col = 0;
			}
			else
			{
				if(tsdf==-1.0)
					col = 100;
				else
					col = 0;
			}
			
			ptr[i]=col;
			//if(tsdf<1.0&&tsdf>-1.0)
			//{
			//	cout<<depth_float[x+y*KINECT_WIDTH]<<endl;
			//	cout<<k<<" "<<j<<" "<<i<<": "<<tsdf<<endl;
			//	//cout<<2.0-i*SIZE_TRIDE<<" "<<j*SIZE_TRIDE+0.5<<": "<<tsdf<<endl;
			//}
		}

		const int nowweight = 1;
		float newtsdf;
		if(tsdfFirstTime)
		{
			newtsdf = tsdf;
			wight = 1;
		}
		else
		{
			unpack_tsdf(pretsdf,preweight,num);
			wight = max(128,preweight+nowweight);
			newtsdf = (pretsdf*preweight+tsdf)/(wight,nowweight);
		}
		pack_tsdf(newtsdf,wight,num);
	}
	//if(k<=g_size/2)
	//{
	//char namebuffer[20];
	//sprintf(namebuffer,"%d.jpg",k);
	//string name = namebuffer;
	//cv::imwrite(name,color);
	//}
	}
	//out.close();
	//cv::pyrUp(color,bigcol,cv::Size(g_size*2,g_size*2));
	//cv::pyrUp(bigcol,bbcol,cv::Size(g_size*4,g_size*4));
	//cv::imshow("zz",bbcol);
	tsdfFirstTime = false;
	cv::imshow("zzz",color);
}
//--------------------------raycast计算--------------------------
float getmaxtime(const ofVec3f &dirc,const ofVec3f &orignal)
{
	float maxx = (dirc.x>0?3.96875-orignal.x:-orignal.x)/dirc.x;
	float maxy = (dirc.y>0?3.96875-orignal.y:-orignal.y)/dirc.y;
	float maxz = (dirc.z>0?3.96875-orignal.z:-orignal.z)/dirc.z;
	return min(min(maxx,maxy),maxz);
}
void testApp::compute_raycast(const ofMatrix4x4 &tMatrix)
{
	ofVec3f LightPosition(0,0,1000);
	cv::Mat rayimage(KINECT_HEIGHT,KINECT_WIDTH,CV_8UC1);

	float tStep = SIZE_TRIDE*0.8;
	ofVec3f startPosition = ofVec3f(-2,-2,0.5);
	int i=0,j=0;
	int num=0;
	int nnum =0;
	for(j=0;j<480;j++)
	{
	uchar *ptr = rayimage.ptr<uchar>(j);
	for(i=0;i<640;i++)
	{
		ptr[i] = 255;
		num = j*KINECT_WIDTH+i;
		pointsmap_final[num] = ofVec3f(0,0,0);
		normalmap_final[num] = ofVec3f(0,0,0);
		ofVec3f rayChange =nmmul(invKcam,ofVec3f(i,j,1));
		///临时变换，以后做T变换
		//rayChange.y *= -1;
		rayChange = tMatrix * rayChange;

		ofVec3f rayOne = 0.5*rayChange - startPosition;
		ofVec3f rayTwo = 1.0*rayChange - startPosition;
		ofVec3f rayFirst = getTheGlobalPostion(rayOne);
		ofVec3f rayEnd = getTheGlobalPostion(rayTwo);
		ofVec3f rayDirct = rayTwo - rayOne;//rayEnd - rayFirst;

		//cout<<rayDirct<<endl;
		rayDirct = rayDirct.normalize(); 
		//cout<<rayDirct<<endl;
		rayDirct.x = (rayDirct.x == 0.f) ? 1e-15 : rayDirct.x;
        rayDirct.y = (rayDirct.y == 0.f) ? 1e-15 : rayDirct.y;
        rayDirct.z = (rayDirct.z == 0.f) ? 1e-15 : rayDirct.z;

		float pretsdf = getTsdfData((int)rayFirst.x,(int)rayFirst.y,(int)rayFirst.z);
		float nowtsdf ;
		float turr = 0;
		float maxtime = getmaxtime(rayDirct,rayOne);
		//cout<<maxtime<<endl;
		ofVec3f rayLast = rayOne + rayDirct*maxtime;
		ofVec3f gLast = getTheGlobalPostion(rayLast);

		//if(getTsdfData((int)gLast.x,(int)gLast.y,(int)gLast.z)<0)
		//{
		ofVec3f rayPre = rayOne;
		while(turr<maxtime)
		{
			ofVec3f rayNext = rayPre + rayDirct*tStep;
			ofVec3f gnum = getTheGlobalPostion(rayNext);
			if(!isvalid(gnum))
			{
				break;
			}
			nowtsdf = getTsdfData((int)gnum.x,(int)gnum.y,(int)gnum.z);
			if(pretsdf>0&&nowtsdf<0)//&&!(pretsdf==1&&nowtsdf==-1))
			{
				//cout<<"first = "<<rayFirst<<endl;
				//cout<<"dircect = "<<rayDirct<<endl;
				//cout<<"raypre = "<<gpre<<"raynext = "<<rayNext<<endl;
				//float sk=pretsdf,skk=nowtsdf;
				//if(pretsdf<1&&nowtsdf>-1)
				//{
				//	cout<<pretsdf<<" "<<nowtsdf<<endl;
				//}

				float Ftdt = triLinary(rayNext);
				if(isnan(Ftdt))
					break;

				float Ft = triLinary(rayPre);
				if(isnan(Ft))
					break;

				float realt = turr - (tStep*Ft/(Ftdt-Ft));
				ofVec3f rayzero = rayPre + rayDirct*(realt-turr);
				pointsmap_final[num] = (rayzero + startPosition);//*100;// 显示用
				//cout<<Ftdt<<" "<<Ft<<endl;
				///计算法向////
				ofVec3f pregnum = getTheGlobalPostion(rayzero);
				if(pregnum.x>1&&pregnum.y>1&&pregnum.z>1&&pregnum.x<TSDF_SIZE-2&&pregnum.y<TSDF_SIZE-2&&pregnum.z<TSDF_SIZE-2)
				{
					ofVec3f t;

					t = rayzero;
					t.x -= SIZE_TRIDE;
					float Fx1 = triLinary(t);
					t = rayzero;
					t.x += SIZE_TRIDE;
					float Fx2 = triLinary(t);
					float nx = (Fx2-Fx1);///(SIZE_TRIDE*2);

					t = rayzero;
					t.y -= SIZE_TRIDE;
					float Fy1 = triLinary(t);
					t = rayzero;
					t.y += SIZE_TRIDE;
					float Fy2 = triLinary(t);
					float ny = (Fy2-Fy1);///(SIZE_TRIDE*2);

					t = rayzero;
					t.z -= SIZE_TRIDE;
					float Fz1 = triLinary(t);
					t = rayzero;
					t.z += SIZE_TRIDE;
					float Fz2 = triLinary(t);
					float nz =(Fz2-Fz1);////(SIZE_TRIDE*2);

					ofVec3f n(nx,ny,nz);
					n = n.normalize();
					normalmap_final[num] = n;
					nnum++;
					
					// 生成图像				
					ofVec3f pp = LightPosition - pointsmap_final[num];
					pp = pp.normalize();
					float r = fabs(pp.dot(n));
					int br = (int)(r*205)+50;
					br = max((int)0,min(br,255));
					ptr[i] = br;
				}
				//num++;
				break;
			}
			else if(pretsdf<0&&nowtsdf>0)
				break;
			else
			{
				pretsdf = nowtsdf;
				rayPre = rayNext;
				turr += tStep;
			}
		}// while()
		//}// if(Last<0)

	}// for(i)
	}// for(j)
	finalPointNum = num;
	//std::cout<<num<<" "<<nnum<<endl;
	cv::imshow("ray",rayimage);

}
//--------------------------PDA与ICP计算------------------------------
ofMatrix3x3 getRotationMatrix(const ofMatrix4x4& Tmaxtrix)
{
	ofMatrix3x3 Rmaxtrix = ofMatrix3x3(Tmaxtrix(0,0),Tmaxtrix(0,1),Tmaxtrix(0,2),
		                               Tmaxtrix(1,0),Tmaxtrix(1,1),Tmaxtrix(1,2),
									   Tmaxtrix(2,0),Tmaxtrix(2,1),Tmaxtrix(2,2));
	//cout<<"Rmaxtrix = "<<endl;
	//cout<<Rmaxtrix<<endl;
	return Rmaxtrix;
}
void testApp::changePosition(const ofMatrix4x4 &changeMatrix,ofVec3f* finalPoints,ofVec3f* finalnorml)
{
	ofMatrix3x3 cnormalMatrix = getRotationMatrix(changeMatrix);

	for(int y = 0; y < KINECT_HEIGHT; y ++) {
		for(int x = 0; x < KINECT_WIDTH; x ++) {
			if(finalPoints[y*KINECT_WIDTH+x].z > 0) {
				ofVec3f point = finalPoints[y*KINECT_WIDTH+x];
				point = changeMatrix * point;
				ofVec3f pnormal = finalnorml[y*KINECT_WIDTH+x];
				pnormal = nmmul(cnormalMatrix,pnormal);
				finalPoints[y*KINECT_WIDTH+x] = point;
				finalnorml[y*KINECT_WIDTH+x] = pnormal;
			}
		}
	}

}
void testApp::compute_icp(float b,const ofVec3f &spacePoint,const ofVec3f &spaceNormal)
{
	Eigen::Vector3f snormal;
	float x = spacePoint.x;
	float y = spacePoint.y;
	float z = spacePoint.z;
	snormal<<spaceNormal.x,spaceNormal.y,spaceNormal.z;
	//cout<<snormal<<endl;
	//cout<<endl;
	Eigen::Matrix<float,3,6> Gx;
	//Gx<<(Eigen::Matrix3f()<<0,-z,y,z,0,-x,-y,x,0).finished(),
	Gx<<(Eigen::Matrix3f()<<0,z,-y,-z,0,x,y,-x,0).finished(),
		 Eigen::MatrixXf::Identity(3,3);
	//cout<<Gx<<endl;
	//cout<<endl;
	Eigen::MatrixXf Gx_t = Gx.transpose();
	Eigen::MatrixXf A_t = Gx_t * snormal;
	Eigen::MatrixXf A = A_t.transpose();
	Eigen::MatrixXf Adata = A_t * A;
	Eigen::MatrixXf Bdata = A_t * b;

	// 相加
	final_A += Adata;
	final_B += Bdata;

}
void testApp::compute_pda(const ofMatrix4x4 &preTmatrix,int timeZ,ofMatrix4x4 &newTmatrix)//,const ofVec3f* pointmapNow,const ofVec3f* pointmap_pre)
{
	ofMatrix4x4 nowTmatrix;
	if(timeZ==0)
		nowTmatrix = preTmatrix;
	else
	{
		// 迭代
		nowTmatrix = newTmatrix;
		//cout<<"z is big than 0ne"<<endl;
	}
	//ofMatrix4x4 invPreT = nowTmatrix.getInverse();
	ofMatrix4x4 invPreT = preTmatrix.getInverse();
	ofMatrix3x3 nowRmatrix = getRotationMatrix(nowTmatrix);
	ofMatrix3x3 preRmatrix = getRotationMatrix(preTmatrix);

	int i=0,j=0;
	int dnum = 0;
	double total_b = 0;
	double total_bb = 0;

	final_A = Eigen::MatrixXf::Zero(6,6);
	final_B = Eigen::MatrixXf::Zero(6,1);

	// 找到对应点对
	for(j=0;j<KINECT_HEIGHT;j++)
	for(i=0;i<KINECT_WIDTH;i++)
	{
		int num = j*KINECT_WIDTH+i;

		if(pointsmap_orignal[num].z>0)//isvalid_orignal[num])
		{
			ofVec4f pointNow = ofVec4f(pointsmap_orignal[num].x,pointsmap_orignal[num].y,pointsmap_orignal[num].z,1);
			ofVec4f spacePointNow = nowTmatrix*pointNow;
			ofVec3f spacenormalNow = nmmul(nowRmatrix,normalmap_orignal[num]);
			//cout<<i<<" "<<j<<endl;

			ofVec4f pointPreC = invPreT*spacePointNow;
			ofVec3f dpointPreC = ofVec3f(pointPreC.x,pointPreC.y,pointPreC.z);
			ofVec3f pointPreZ = nmmul(the_K_cam,dpointPreC);
			ofVec2f pointPreU = ofVec2f(pointPreZ.x/pointPreZ.z,pointPreZ.y/pointPreZ.z);
			//cout<<endl;

			int num_p = (int)(pointPreU.x+0.5)+(int)(pointPreU.y+0.5)*KINECT_WIDTH;
			int kk = num;
			if(pointsmap_final[num_p].z>0)
			{
				ofVec4f pointPre = ofVec4f(pointsmap_final[num_p].x,pointsmap_final[num_p].y,pointsmap_final[num_p].z,1);
				ofVec4f spacePointPre = pointPre;
				//ofVec4f spacePointPre = preTmatrix*pointPre;
				ofVec3f spacenormalPre = normalmap_final[num_p];
				//ofVec3f spacenormalPre = nmmul(preRmatrix,normalmap_final[num_p]);
				float pointDistance = spacePointNow.distance(spacePointPre);
				float pointAngle = spacenormalNow.angle(spacenormalPre);
				if(pointDistance<=THRESHOLD_D&&pointAngle<=THRESHOLD_A)
				{
					//cout<<spacenormalPre<<" "<<spacenormalNow<<endl;
					// 加入点对
					//if(pointDistance>0.001)
					//{
					//	cout<<i<<" "<<j<<endl;
					//	cout<<spacePointNow<<endl;
					//	cout<<spacePointPre<<endl;
					//}
					//spacePointNow = newTmatrix * spacePointNow;
					ofVec3f dspacePointNow = ofVec3f(spacePointNow.x,spacePointNow.y,spacePointNow.z);
					ofVec3f dspacePointPre = ofVec3f(spacePointPre.x,spacePointPre.y,spacePointPre.z);
					ofVec3f minsPoint = (dspacePointPre - dspacePointNow);
					float b =  minsPoint.dot(spacenormalPre);
					if(b > 0)
						total_b += b;
					else
						total_bb += b;
					compute_icp(b,dspacePointNow,spacenormalPre);
					++dnum;
				}
			}
		}
	}
	cout<<"dnum = "<<dnum<<" total_b = "<<total_b<<" total_bb = "<<total_bb<<endl;
	//cout<<"final a ="<<endl;
	//cout<<final_A<<endl;
	//cout<<"final b ="<<endl;
	//cout<<final_B<<endl;

	//double deta = final_A.determinant();
	//if(deta>0) // 条件待修改
		Eigen::Matrix<float,6,1> result = final_A.llt().solve(final_B).cast<float>();
		
	cout<<"result = "<<endl;
	cout<<result<<endl;

	float alpha = result(0);
	float beta = result(1);
	float gamma = result(2);

	Eigen::Matrix3f rinc = (Eigen::Matrix3f)Eigen::AngleAxisf(gamma,Eigen::Vector3f::UnitZ())*Eigen::AngleAxisf(beta,Eigen::Vector3f::UnitY())*Eigen::AngleAxisf(alpha,Eigen::Vector3f::UnitX());
	//Eigen::Matrix3f rinc = (Eigen::Matrix3f)Eigen::AngleAxisf(alpha,Eigen::Vector3f::UnitX())*Eigen::AngleAxisf(beta,Eigen::Vector3f::UnitY())*Eigen::AngleAxisf(gamma,Eigen::Vector3f::UnitZ());
	//cout<<"rinc = "<<endl;
	//cout<<rinc<<endl;

	ofMatrix4x4 changeMatrix = ofMatrix4x4(rinc(0,0),rinc(0,1),rinc(0,2),result(3),
		                                   rinc(1,0),rinc(1,1),rinc(1,2),result(4),
		                                   rinc(2,0),rinc(2,1),rinc(2,2),result(5),
										   0,0,0,1);
	cout<<"change = "<<endl;
	cout<<changeMatrix<<endl;

	newTmatrix = changeMatrix * newTmatrix;
	cout<<"newTM = "<<endl;
	cout<<newTmatrix<<endl;

	//changePosition(newTmatrix,pointsmap_orignal,normalmap_orignal);
	//changePosition(newTmatrix,pointsmap_final);
}
//---------------------------------------------------------------
void testApp::update() {
	ofBackground(100, 100, 100);
	NUI_IMAGE_FRAME pImageFrame = {0};
	NUI_LOCKED_RECT locked_rect = {0};
	if(bKinectInitSucessful)
	{
		HRESULT hr=m_nui->NuiImageStreamGetNextFrame(depthStreamHandle,0,&pImageFrame);
		if (SUCCEEDED(hr))
		{
			hr=pImageFrame.pFrameTexture->LockRect(0,&locked_rect,NULL,0);
			if(locked_rect.Pitch!=0) {
		
			// load grayscale depth image from the kinect source
			//grayImage.setFromPixels(locked_rect.pBits, KINECT_WIDTH, KINECT_HEIGHT);
			// 将数据存入Mat中
			cv::Mat depthimage(KINECT_HEIGHT,KINECT_WIDTH,CV_8UC1);
			//cv::Mat depthimage2(KINECT_HEIGHT,KINECT_WIDTH,CV_8UC1);

			if(test)
			{
			for (int i=0; i<depthimage.rows; i++) 
			{
				uchar *ptr = depthimage.ptr<uchar>(i);  //第i行的指针
				//uchar *ptr2 = depthimage2.ptr<uchar>(i);
				uchar *pBufferRun = (uchar*)(locked_rect.pBits) + i * locked_rect.Pitch;
				USHORT *pBuffer = (USHORT*) pBufferRun;					 
				for (int j=0; j<depthimage.cols; j++) 
				{
					//kinect sdk 方法求三维点坐标
					//if(pBuffer[j]>0)
					//{
					//	//cout<<i<< " "<<j<<" "<<vec.x*1000<<","<<vec.y*1000<<","<<vec.z*1000<<endl;
					//	Vector4 vec = NuiTransformDepthImageToSkeleton(j,i,pBuffer[j],NUI_IMAGE_RESOLUTION_640x480);
					//	ofVec3f of(vec.x*500,vec.y*500,vec.z*500);
					//	pointsmap_orignal[i*depthimage.cols+j] = of;
					//}
					//else
					//{
					//	pointsmap_orignal[i*depthimage.cols+j] = ofVec3f(0.0,0.0,0.0);
					//}

					pBuffer[j]=pBuffer[j]>>3;
					//if(pBuffer[j]>0&&j<10)
					//{
					//	cout<<(pBuffer[j])<<" ";
					//}
					//depth_float[i*depthimage.cols+j]=(float)(pBuffer[j])/10;
					ptr[j]=(uchar)(256*pBuffer[j]/0x0fff); //直接将数据归一化处理
					//cout<<depth_float[i*depthimage.cols+j]<<" ";
					//ptr2[j]=(uchar)(256*depth_float[i*depthimage.cols+j]/0x0fff);
				} 
			}
			//test=false;
			}
			// Mat类型转换 耗时（废弃）
			//depthimage.convertTo(depthimage,CV_32F);

			const int photonumMax = 1;
			// 滤波 双边滤波效果没出现（？）
			cv::imshow("orginal",depthimage);
			//cv::imshow("sss",depthimage2);			
			cv::Mat depthimage_filter = depthimage.clone();
			// 读取照片
#if 0
			if(photonum < photonumMax)
			{
				char namebuffer[20];
				sprintf(namebuffer,"%d.jpg",photonum);
				string name = namebuffer;
				depthimage_filter = cv::imread(name,CV_LOAD_IMAGE_GRAYSCALE);
				++photonum;
			}
#endif
			//cv::medianBlur(depthimage,depthimage_filter,5);
			//clock_t t1=clock();
			cv::bilateralFilter(depthimage,depthimage_filter,11,20,20);
			//clock_t t2=clock();
			//double ti=(double)(t2-t1)/CLOCKS_PER_SEC;
			//cout<<ti<<endl;
			//cv::GaussianBlur(depthimage,depthimage_filter,cv::Size(5,5),0,0);
			//cv::imshow("filter",depthimage_filter);
			// 保存照片
#if 0 
			if(photonum < photonumMax)
			{
				char namebuffer[20];
				sprintf(namebuffer,"%d.jpg",photonum);
				string name = namebuffer;
				cv::imwrite(name,depthimage);
				++photonum;
			}
#endif
			// pyramid creat 高斯降采 论文中提到避免边界被平滑未处理（？）
		
			//cv::Mat downmap1;
			//cv::pyrDown(depthimage_filter,downmap1,cv::Size(depthimage_filter.cols/2,depthimage_filter.rows/2));
			//cv::imshow("downonce",downmap1);
			//cv::Mat downmap2;
			//cv::pyrDown(downmap1,downmap2,cv::Size(downmap1.cols/2,downmap1.rows/2));
			//cv::imshow("downtwice",downmap2);
		
			// 滤波后深度数据提取与单位转换
			
			//ofstream out;
			//if(photonum < photonumMax&&test3)
			//{
			//	char namebuffer[10];
			//	sprintf(namebuffer,"%d.txt",photonum);
			//	out.open(namebuffer,ios::trunc);
			//}
			for (int i=0; i<depthimage_filter.rows; i++) 
			{
				uchar *ptr = depthimage_filter.ptr<uchar>(i);  //第i行的指针			 
				for (int j=0; j<depthimage_filter.cols; j++) 
				{
					float s = (float)(ptr[j]*0x0fff)/256;
					//if(s>2500&&s<2600)
					//{
					//	cout<<s<<endl;
					//}
					depth_float[i*depthimage_filter.cols+j] = s/1000;//显示点云用，计算时要除1000，化为米单位
					//if(test3)
					//	out<<s<<" ";				
				}
			}
			//out.close();
			test3 = false;
			// 试验区
			ofMatrix4x4 tk=ofMatrix4x4(1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1);
			ofMatrix4x4 tk_next = tk;
			ofVec4f camp=ofVec4f(0,0,0,1);
			int size_g=TSDF_SIZE;
			const int maxcountnum = 3;
			if(countnum < maxcountnum)
			{
				//compute_tsdf(tk,depthimage,size_g,camp);
				++countnum;
				if(countnum == maxcountnum)
					test2 = true;
			//}
			//if(test2)
			//{
				//cv::imwrite("original.jpg",depthimage_filter);
				//cout<<tk<<endl;
				clock_t t1=clock();
				compute_points(depthimage_filter,depthimage_filter.rows,depthimage_filter.cols,pointsmap_orignal);
				compute_tsdf(tk_next,depthimage,size_g,camp);
				compute_raycast(tk_next);
				compute_normal(pointsmap_orignal,depthimage_filter.rows,depthimage_filter.cols,normalmap_orignal);
				//compute_normal(pointsmap_final,depthimage_filter.rows,depthimage_filter.cols,normalmap_final);

				for(int i = 0;i < 10; ++ i)
				{
					compute_pda(tk,i,tk_next);
				}
				camp = tk_next * camp;
				tk = tk_next;
				//changePosition(tk_next,pointsmap_orignal,normalmap_orignal);
				clock_t t2=clock();
				double ti=(double)(t2-t1)/CLOCKS_PER_SEC;
				std::cout<<ti<<endl;
				//test2=false;
			}

			// 转换坐标系为摄像机坐标系与法向

			//changeType(depthimage_filter,depthimage_filter.rows,depthimage_filter.cols);//无用了
			//compute_points(depthimage_filter,depthimage_filter.rows,depthimage_filter.cols,pointsmap_orignal);
			//compute_points(downmap1,downmap1.rows,downmap1.cols,pointsmap_downonce);
			//compute_points(downmap2,downmap2.rows,downmap2.cols,pointsmap_downtwice);
			//
			//compute_normal(pointsmap_orignal,depthimage_filter.rows,depthimage_filter.cols,normalmap_orignal);
			//compute_normal(pointsmap_downonce,downmap1.rows,downmap1.cols,normalmap_downonce);
			//compute_normal(pointsmap_downtwice,downmap2.rows,downmap2.cols,normalmap_downtwice);

			// 其他计算矩阵乘法
			// 方法1
			//float array1[]={1,0,1,0,2,1,0,0,3};
			//float array2[]={1,0,1};
			//CvMat *A=cvCreateMat(3,3,CV_32FC1);
			//CvMat *B=cvCreateMat(3,1,CV_32FC1);
			//CvMat *ResultMatrix=cvCreateMat(3,1,CV_32FC1);
			//cvSetData(A,array1,A->step);
			//cvSetData(B,array2,B->step);
			//cvMatMul(A,B,ResultMatrix);
			//cv::Mat ssa=cv::Mat(ResultMatrix,true);
			//cout<<ssa<<endl;
			// 方法2
			//cv::Mat cam_k=(cv::Mat_<double>(3,3)<<1,0,0,0,2,0,0,0,3);
			//cv::Mat L=(cv::Mat_<double>(3,1)<<1,1,0);
			//cv::Mat result=cv::Mat(3,1,CV_32F);
			//cout<<L<<endl;
			//cvmMul(cam_k,L,result);
			//cvMatMul(&cam_k,&L,&result);
			//CvMat a=cam_k;
			//CvMat b=L;
			//CvMat c=result;
			//cvMatMul(a,b,c);
		
			// we do two thresholds - one for the far plane and one for the near plane
			// we then do a cvAnd to get the pixels which are a union of the two thresholds
			if(bThreshWithOpenCV) { 
				//grayThreshNear = grayImage;
				//grayThreshFar = grayImage;
				//grayThreshNear.threshold(nearThreshold, true);
				//grayThreshFar.threshold(farThreshold);
				//cvAnd(grayThreshNear.getCvImage(), grayThreshFar.getCvImage(), grayImage.getCvImage(), NULL);
			} else {
				//grayImage.setFromPixels(depthimage_filter.data, KINECT_WIDTH, KINECT_HEIGHT);
				// or we do it ourselves - show people how they can work with the pixels
				//int numPixels = grayImage.getWidth() * grayImage.getHeight();
				//for(int i = 0; i < numPixels; i++) {
				//	if(pix[i] < nearThreshold && pix[i] > farThreshold) {
				//		pix[i] = 255;
				//	} else {
				//		pix[i] = 0;
				//	}
				//}
			}
			}
			m_nui->NuiImageStreamReleaseFrame(depthStreamHandle,&pImageFrame);
		}
	}
}

#pragma region draw_function
//--------------------------------------------------------------
void testApp::draw() {
	
	ofSetColor(255, 255, 255);
	
	if(bDrawPointCloud) {
		easyCam.begin();
		drawPointCloud();
		easyCam.end();
	} else {		
		//grayImage.draw(10, 320, 400, 300);
		//contourFinder.draw(10, 320, 400, 300);
	}
	
	// draw instructions
	ofSetColor(255, 255, 255);
	stringstream reportStream;
           
	reportStream << "press p to switch between images and point cloud, rotate the point cloud with the mouse" << endl
	<< "using opencv threshold = " << bThreshWithOpenCV <<" (press spacebar)" << endl
	<< "set near threshold " << nearThreshold << " (press: + -)" << endl
	<< "set far threshold " << farThreshold << " (press: < >) num blobs found " 
	<< ", fps: " << ofGetFrameRate() << endl;
    
	ofDrawBitmapString(reportStream.str(), 20, 652);    
}

void testApp::drawPointCloud() {
	int w = 640;
	int h = 480;
	float x_avg=0;
	float y_avg=0;
	float z_avg=0;
	int valid_num=0;
	//ofMesh mesh;
	mesh.setMode(OF_PRIMITIVE_POINTS);
	//mesh.setMode(OF_PRIMITIVE_TRIANGLES);
	int step = 4;
	ofColor balck;
	balck.r=0;
	balck.g=0;
	balck.b=0;

	ofColor white;
	white.r=255;
	white.g=255;
	white.b=255;

	if(test)
	{
		mesh.clear();
		if(bThreshWithOpenCV)
		{
			//for(int y = 0; y < h; y+= step)
			//{
			//	mesh.addColor(balck);
			//	ofVec3f of(y,0,1200);
			//	mesh.addVertex(of);
			//}

			for(int y = 0; y < h; y += step) {
				for(int x = 0; x < w; x += step) {
					if(pointsmap_orignal[y*w+x].z > 0) {
						mesh.addColor(balck);
						//pointsmap_orignal[y*w+x].y *= -1;
						mesh.addVertex(pointsmap_orignal[y*w+x]*100);
						mesh.addNormal(normalmap_orignal[y*w+x]);
					}
				}
			}

			for(int y = 0; y < h; y += step) {
				for(int x = 0; x < w; x += step) {
					if(pointsmap_final[y*w+x].z > 0) {
						mesh.addColor(balck);
						//pointsmap_final[y*w+x].y *= -1;
						mesh.addVertex(pointsmap_final[y*w+x]*100);
						mesh.addNormal(normalmap_final[y*w+x]);
					}
				}
			}
			//for(int mm = 0;mm<finalPointNum;mm += 1)
			//{
			//	mesh.addColor(balck);
			//	mesh.addVertex(pointsmap_final[mm]);
			//	mesh.addNormal(normalmap_final[mm]);
			//}
		}
		else
		{
			//int row=0;
			//int col=0;
			//int tt=0;
			//unsigned char * pix = grayImage.getPixels();
			//for(int y = 0; y < h; y += step) {
			//	for(int x = 0; x < w; x += step) {			
			//		row=y/2;
			//		col=x/2;
			//		ofVec3f a(0,0,0);
			//		if(pix[x+w*y]>0)//farThreshold)
			//		{
			//			//ss[row][col]=1;
			//			////mesh.addColor(balck);
			//			//mesh.addColor(kinect.getColorAt(x,y));
			//			//a=kinect.getWorldCoordinateAt(x,y);
			//			//valid_num++;
			//			//x_avg+=a.x;
			//			//y_avg+=a.y;
			//			//z_avg+=a.z;
			//			////cout<<a.x<<" "<<a.y<<" "<<a.z<<endl;
			//			//mesh.addVertex(kinect.getWorldCoordinateAt(x,y));
			//		}
			//		else
			//		{
			//			ss[row][col]=0;
			//			//mesh.addVertex(a);
			//			//mesh.addColor(white);
			//		}
			//		//my_points.push_back(a);
			//		//mesh.addVertex(a);
			//		//m_points[col+row*320]=a;
			//		//m_points[x+y*320]=a;
			//		//mesh.addVertex(kinect.getWorldCoordinateAt(x,y));			
			//	}
			//}
			////cout<<mesh.getNumVertices()<<endl;
			//x_avg/=valid_num;
			//y_avg/=valid_num;
			//z_avg/=valid_num;
			////cout<<"avg:"<<x_avg<<" "<<y_avg<<" "<<z_avg<<endl;
		}
		test=false;
	}

	glPointSize(3);
	ofPushMatrix();
	// the projected points are 'upside down' and 'backwards' 
	ofEnableDepthTest();
	//x轴为红，y轴为绿，z轴为蓝
	ofDrawAxis(100);

	ofScale(1, 1, -1);
	ofTranslate(x_back, y_back,z_back); // center the points a bit
	//ofEnableDepthTest();
	mesh.drawVertices();
	if(draww)
	{
		vector<ofVec3f> point_map = mesh.getVertices();
		vector<ofVec3f> normal_map = mesh.getNormals();
		ofSetColor(ofColor::red);
		ofSetLineWidth(1);
		vector<ofVec3f>::iterator it;
		vector<ofVec3f>::iterator is;
		for(it = point_map.begin(),is=normal_map.begin();it < point_map.end()&&is< normal_map.end();it++,is++)
		{
			ofLine(it->x,it->y,it->z,it->x+3*is->x,it->y+3*is->y,it->z+3*is->z);
		}

		//for(int y = 0; y < h; y += step) {
		//		for(int x = 0; x < w; x += step) {
		//			if(pointsmap_orignal[y*w+x].z > 0) {
		//				ofVec3f point_a=pointsmap_orignal[y*w+x];
		//				ofVec3f normal_a=normalmap_orignal[y*w+x];
		//				ofLine(point_a.x,point_a.y,point_a.z,point_a.x+3*normal_a.x,point_a.y+3*normal_a.y,point_a.z+3*normal_a.z);
		//			}
		//		}
		//	}
	
	}

	//if(bThreshWithOpenCV)
	//	mesh.drawVertices();
	//else
	//{
	//	ofVec3f* points_line=m_points;
	//	ofVec3f* points_next_line=m_points+320;
	//	bool mesh_break=true;
	//	for(int y = 0; y < 240-1; y++) {
	//		for(int x = 0; x < 320; x++) {
	//			ofVec3f& s_point1=points_line[x];
	//			ofVec3f& s_point2=points_next_line[x];
	//			if(((s_point1.z<400||s_point1.z>1000)||(s_point2.z<400||s_point2.z>1000)))
	//				//||((s_point1.x<x_avg-200||s_point1.x>x_avg+200)||(s_point2.x<x_avg-200||s_point2.x>x_avg+200))||((s_point1.y<y_avg-200||s_point1.y>y_avg+200)||(s_point2.y<y_avg-200||s_point2.y>y_avg+200)))
	//			{
	//				if(!mesh_break)
	//				{
	//					mesh_break=true;
	//					glEnd();
	//				}
	//				continue;
	//			}
	//			if(mesh_break)
	//			{
	//				if(draww)
	//				{
	//					glBegin(GL_TRIANGLE_STRIP);
	//				}
	//				else
	//				{
	//				   glBegin(GL_POINTS);
	//				}
	//				mesh_break=false;
	//			}
	//			
	//			glColor3f(1.0,1.0,1.0);
	//			glVertex3f(s_point1.x,s_point1.y,s_point1.z);
	//			glColor3f(1.0,1.0,1.0);
	//			glVertex3f(s_point2.x,s_point2.y,s_point2.z);
	//		}
	//		if(!mesh_break)
	//		{
	//			glEnd();
	//			mesh_break=true;
	//		}
	//		points_line+=320;
	//		points_next_line+=320;
	//	}
	//}
	//mesh.draw();
	ofDisableDepthTest();
	ofPopMatrix();
}
#pragma endregion
//--------------------------------------------------------------
void testApp::exit() {

}

//--------------------------------------------------------------
void testApp::keyPressed (int key) {
	switch (key) {
		case ' ':
			bThreshWithOpenCV = !bThreshWithOpenCV;
			test=true;
			break;
			
		case'p':
			bDrawPointCloud = !bDrawPointCloud;
			break;
			
		case '>':
		case '.':
			farThreshold ++;
			if (farThreshold > 255) farThreshold = 255;
			break;
			
		case '<':
		case ',':
			farThreshold --;
			if (farThreshold < 0) farThreshold = 0;
			break;
			
		case '+':
		case '=':
			nearThreshold ++;
			if (nearThreshold > 255) nearThreshold = 255;
			break;
			
		case '-':
			nearThreshold --;
			if (nearThreshold < 0) nearThreshold = 0;
			break;
						
		case 'n':
			draww=!draww;
			break;

		case 'w':
			++ y_back;
		break;

		case 's':
			-- y_back;
		break;

		case 'a':
			-- x_back;
		break;

		case 'd':
			++ x_back;
		break;

		case 'x':
			-- z_back;
		break;

		case 'z':
			++ z_back;
		break;
	}
}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button)
{}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button)
{}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button)
{}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h)
{}
