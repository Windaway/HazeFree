#ifndef HazeREMOVE_H
#define HazeREMOVE_H
//算法需要一个uchar行向量的图像做为输入，RGB排列，具体参见例子使用
int AHR(unsigned char* pImage, unsigned int uiXRes, unsigned int uiYRes, unsigned char Min,
	unsigned char Max, unsigned int uiNrX, unsigned int uiNrY,
	unsigned int uiNrBins, float fCliplimit);
	/*	融合去雾算法
	*	pImage - 输入图像指针
	*   uiXRes - 宽
	*   uiYRes - 高
	*   Min - 最低灰度级0
	*   Max - 最高灰度级255
	*   uiNrX - 横向分区数 算法不做padding
	*   uiNrY - 纵向分区数 算法不做padding，需要设置为可以整除
	*   uiNrBins - 直方图数量/256默认，当然如果弄128速度更快,而且效果也很好的。
	*   fCliplimit - 规则化的截断下限，1则为无对比度限制的自适应直方图均衡，一般设成0.01即可
	*/
int CLAHE1(unsigned char* pImage,  unsigned int uiXRes, unsigned int uiYRes, unsigned char Min,
	unsigned char Max, unsigned int uiNrX, unsigned int uiNrY,
	unsigned int uiNrBins, float fCliplimit);
	/*CLAHE去雾算法*/
int DPHR(unsigned char* pImage, unsigned int uiXRes, unsigned int uiYRes);
	/*暗通道去雾算法*/

#endif