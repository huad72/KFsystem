#pragma once

#include "ofMain.h"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <Windows.h>
#include <Ole2.h>
#include "NuiApi.h"

#include <Eigen\Core>
#include <Eigen\SVD>
#include <Eigen\Cholesky>
#include <Eigen\Geometry>
#include <Eigen\LU>
//#include <algorithm>
// uncomment this to read from two kinects simultaneously
//#define USE_TWO_KINECTS

struct short2{
	short tsdf;
	short weight;
};

class testApp : public ofBaseApp {
public:
	
	void setup();
	void update();
	void draw();
	void exit();
	
	void drawPointCloud();
	
	/////myadd/////
	void testh();
	void readglobalfile(char *filename,float *data);
	// 基础计算
	ofVec3f nmmul(const ofMatrix3x3 &mat1,const ofVec3f &mat2) const;
	void changeType(const cv::Mat &depthimage,const int rows,const int cols);
	void compute_normal(const ofVec3f *vec,const int rows,const int cols,ofVec3f *normalmap);
	void compute_points(const cv::Mat &depthimage,const int rows,const int cols,ofVec3f *pointsmap);

	// tsdf计算
	void pack_tsdf(const float &tsdf,const int &weight,const int num);
	void unpack_tsdf(float &tsdf,int &weight,const int num);
	void compute_tsdf(const ofMatrix4x4 &T_G,cv::Mat &depthimage,const int g_size,const ofVec4f &camP);
	// raycast计算
	ofVec3f getTheGlobalPostion(const ofVec3f &gPoint);
	float getTsdfData(const int x,const int y,const int z);
	float triLinary(const ofVec3f &point);
	void compute_raycast();
	// ICP计算
	void compute_pda(const ofMatrix4x4 &preTmatrix,int timeZ,ofMatrix4x4 &newTmatrix);
	void compute_icp(float b,const ofVec3f &spacePoint,const ofVec3f &spaceNormal);
	void changePosition(const ofMatrix4x4 &changeMatrix,ofVec3f* finalPoints,ofVec3f* finalnorml);


	void keyPressed(int key);
	void mouseDragged(int x, int y, int button);
	void mousePressed(int x, int y, int button);
	void mouseReleased(int x, int y, int button);
	void windowResized(int w, int h);
	
	//ofxKinect kinect;
	HANDLE nextColorFrameEvent;
	HANDLE depthStreamHandle;
	INuiSensor *m_nui;
	
	//ofxCvColorImage colorImg;
	//
	//ofxCvGrayscaleImage grayImage; // grayscale depth image
	//ofxCvGrayscaleImage grayThreshNear; // the near thresholded image
	//ofxCvGrayscaleImage grayThreshFar; // the far thresholded image
	//
	//ofxCvContourFinder contourFinder;
	bool bKinectInitSucessful;
	bool bThreshWithOpenCV;
	bool bDrawPointCloud;
	
	int nearThreshold;
	int farThreshold;
	
	int angle;
	
	// used for viewing the point cloud
	ofEasyCam easyCam;
	ofMesh mesh;

	// 我添加的变量
	bool test;
	bool test2;
	bool test3;
	bool draww;
	int countnum;// tsdf多帧读取控制
	int photonum;// 拍摄照片
	int ss[245][325];
	int finalPointNum;

	ofVec3f m_points[325000];

	bool tsdfFirstTime;
	float *depth_float;
	float *tsdfData;
	short2 *tsdfWeightData;

	Eigen::MatrixXf final_A;
	Eigen::MatrixXf final_B;

	ofMatrix3x3 the_K_cam;
	ofMatrix3x3 invKcam;
	ofVec3f normalmap_orignal[325000];
	ofVec3f normalmap_downonce[90000];
	ofVec3f normalmap_downtwice[20000];

	ofVec3f pointsmap_orignal[325000];
	ofVec3f pointsmap_downonce[90000];
	ofVec3f pointsmap_downtwice[20000];

	ofVec3f pointsmap_final[325000];
	ofVec3f normalmap_final[325000];

	bool isvalid_orignal[325000];
	//bool isvalid_downonce[90000];
	//bool isvalid_downtwice[20000];

	//vector<ofVec3f> my_points;
};
inline ofVec3f testApp::nmmul(const ofMatrix3x3 &mat1,const ofVec3f &mat2) const{

	ofVec3f resultmat(0);
	resultmat.x=mat1.a*mat2.x+mat1.b*mat2.y+mat1.c*mat2.z;
	resultmat.y=mat1.d*mat2.x+mat1.e*mat2.y+mat1.f*mat2.z;
	resultmat.z=mat1.g*mat2.x+mat1.h*mat2.y+mat1.i*mat2.z;
	return resultmat;
}