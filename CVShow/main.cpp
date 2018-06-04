#include "stdio.h"
#include <opencv2/opencv.hpp>
#include "hazeremove.h"
using namespace cv;

int main(void){

	char* imagename = "forest.jpg"; //输入文件名并将图像写入内存
	
	Mat img = imread(imagename);
	Mat dst;
	resize(img, dst, Size(1280, 720), (0, 0), (0, 0), INTER_LINEAR);
	uchar* mPtr = new uchar[dst.total()*3];
	int a = dst.total();
	memcpy(mPtr, dst.data, dst.total() * 3);

	
	AHR(mPtr, 1280, 720, 0, 255, 8, 8, 256, 0.01);
	//CLAHE1(mPtr, 1280, 720, 0, 255, 8, 8, 256, 0.01);
	//DPHR(mPtr, 1280, 720);
	//可以分别看三种结果的效果图
	memcpy(dst.data, mPtr , dst.total() * 3);
	imshow("haha",dst);
	waitKey(0);
	return 0;
}