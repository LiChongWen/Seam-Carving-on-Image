
#include "sc.h"
#include "math.h"
#include <algorithm>
#include <vector>
#include <string>

using namespace cv;
using namespace std;

class Vertex{
public:
    int x;
    int y;
    float energy;

    Vertex(){}
    Vertex(int x, int y, float energy){
        this->x = x;
        this->y = y;
        this->energy = energy;
    }
    int connectedX;
    int connectedY;

};

float cGray(int r, int g, int b) {
    return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

void calculateEnergry(Mat &gray, Mat &energy, int rows, int cols){

    for (int i = 1; i < rows-1; ++i) {
        for (int j = 1; j < cols-1; ++j) {
            float x = gray.at<float>(i, j-1) - gray.at<float>(i, j+1);
            float y = gray.at<float>(i-1, j) - gray.at<float>(i+1, j);
            energy.at<float>(i, j) = sqrt( (x*x) + (y*y) );
            if(i==1){
                float fx = gray.at<float>(0, j-1) - gray.at<float>(0, j+1);//第一行
                float fy = gray.at<float>(1, j);
                energy.at<float>(0,j) = sqrt((fx*fx)+(fy*fy));

                float lx = gray.at<float>(rows-1, j-1) - gray.at<float>(rows-1, j+1);//last一行
                float ly = gray.at<float>(rows-2, j);
                energy.at<float>(rows-1,j) = sqrt((lx*lx)+(ly*ly));
            }
        }
        float fx = gray.at<float>(i, 1);
        float fy = gray.at<float>(i-1, 0) - gray.at<float>(i+1, 0);
        energy.at<float>(i, 0) = sqrt( (fx*fx) + (fy*fy) ); //第一列

        float lx = gray.at<float>(i, cols-2);
        float ly = gray.at<float>(i-1, cols-1) - gray.at<float>(i+1, cols-1);
        energy.at<float>(i, cols-1) = sqrt( (lx*lx) + (ly*ly) ); //last一列
    }
    energy.at<float>(0, 0) = sqrt(energy.at<float>(0, 1) * energy.at<float>(1, 0));
    energy.at<float>(rows-1, cols-1) = sqrt(energy.at<float>(rows-1, cols-2) * energy.at<float>(rows-2, cols-1));
    energy.at<float>(rows-1, 0) = sqrt(energy.at<float>(rows-2 , 0)* energy.at<float>(rows-1, 1));
    energy.at<float>(0, cols-1) = sqrt(energy.at<float>(0, cols-2)*energy.at<float>(1, cols-1));
}

void optimalVEnergy(Mat &energy, vector<vector<Vertex > > &oEnergy, int rows, int cols){
    //vector<vector<Vertex > > oEnergy = oEnergy1;
    Vertex empty = Vertex(-1,-1,-1);
    for (int l = 0; l < rows; ++l) {
        for (int i = 0; i < cols; ++i) {
            oEnergy[l][i].x = l;
            oEnergy[l][i].y = i;
        }
    }
    for (int j = 0; j < cols; ++j) {
        oEnergy[0][j].energy = energy.at<float>(0, j);
        oEnergy[0][j].connectedX = -1;
        oEnergy[0][j].connectedY = -1;
        if(oEnergy[0][j].energy<0){
            int y = j;
            bool pass = true;
        }
    }
    for (int i = 1; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if(j==0 || j == cols-1){
                if(j==0){
                    if(oEnergy[i-1][j].energy < oEnergy[i-1][j+1].energy){
                        oEnergy[i][j].energy = energy.at<float>(i, j) + oEnergy[i-1][j].energy;
                        oEnergy[i][j].connectedX = i-1;
                        oEnergy[i][j].connectedY = j;
                        if(oEnergy[i][j].energy<0){
                            int x = i;
                            int y = j;
                            bool pass = true;
                        }
                    } else{
                        oEnergy[i][j].energy = energy.at<float>(i, j) + oEnergy[i-1][j+1].energy;
                        oEnergy[i][j].connectedX = i-1;
                        oEnergy[i][j].connectedY = j+1;
                        if(oEnergy[i][j].energy<0){
                            int x = i;
                            int y = j;
                            bool pass = true;
                        }
                    }
                } else{
                    if(oEnergy[i-1][j].energy < oEnergy[i-1][j-1].energy){
                        oEnergy[i][j].energy = energy.at<float>(i, j) + oEnergy[i-1][j].energy;
                        oEnergy[i][j].connectedX = i-1;
                        oEnergy[i][j].connectedY = j;
                        if(oEnergy[i][j].energy<0){
                            int x = i;
                            int y = j;
                            bool pass = true;
                        }
                    } else{
                        oEnergy[i][j].energy = energy.at<float>(i, j) + oEnergy[i-1][j-1].energy;
                        oEnergy[i][j].connectedX = i-1;
                        oEnergy[i][j].connectedY = j-1;
                        if(oEnergy[i][j].energy<0){
                            int x = i;
                            int y = j;
                            bool pass = true;
                        }
                    }
                }
            }else{

                float minV = min(oEnergy[i-1][j-1].energy, min(oEnergy[i-1][j].energy, oEnergy[i-1][j+1].energy));
                for (int k = j-1; k <= j+1; ++k) {
                    if(oEnergy[i-1][k].energy == minV){
                        oEnergy[i][j].energy = energy.at<float>(i, j) + minV;
                        oEnergy[i][j].connectedX = i-1;
                        oEnergy[i][j].connectedY = k;
                        if(oEnergy[i][j].energy<0){
                            int x = oEnergy[i-1][j-1].energy;
                            int y = oEnergy[i-1][j].energy;
                            int z = oEnergy[i][j+1].energy;
                            bool pass = true;
                        }
                        break;
                    }
                }
            }
        }
    }

    //oEnergy1 = oEnergy;
}

void optimalHEnergy(Mat &energy, vector<vector<Vertex > > &oEnergy, int rows, int cols){
    Vertex empty = Vertex(-1,-1,-1);
    for (int l = 0; l < rows; ++l) {
        for (int i = 0; i < cols; ++i) {
            oEnergy[l][i].x = l;
            oEnergy[l][i].y = i;
        }
    }
    for (int i = 0; i < rows; ++i) {
        oEnergy[i][0].energy = energy.at<float>(i, 0);
        oEnergy[i][0].connectedX = -1;
        oEnergy[i][0].connectedY = -1;
    }
    for (int j = 1; j < cols; ++j) {
        for (int i = 0; i < rows; ++i) {
            if(i==0 || i == rows-1){
                if(i==0){
                    if(oEnergy[i][j-1].energy < oEnergy[i+1][j-1].energy){
                        oEnergy[i][j].energy = energy.at<float>(i, j) + oEnergy[i][j-1].energy;

                        oEnergy[i][j].connectedX = i;
                        oEnergy[i][j].connectedY = j-1;
                    } else{
                        oEnergy[i][j].energy = energy.at<float>(i, j) + oEnergy[i+1][j-1].energy;

                        oEnergy[i][j].connectedX = i+1;
                        oEnergy[i][j].connectedY = j-1;
                    }
                } else{
                    if(oEnergy[i][j-1].energy < oEnergy[i-1][j-1].energy){
                        oEnergy[i][j].energy = energy.at<float>(i, j) + oEnergy[i][j-1].energy;

                        oEnergy[i][j].connectedX = i;
                        oEnergy[i][j].connectedY = j-1;
                    } else{
                        oEnergy[i][j].energy = energy.at<float>(i, j) + oEnergy[i-1][j-1].energy;

                        oEnergy[i][j].connectedX = i-1;
                        oEnergy[i][j].connectedY = j-1;
                    }
                }
            }else{
                float minV = min(oEnergy[i-1][j-1].energy, min(oEnergy[i][j-1].energy, oEnergy[i+1][j-1].energy));
                for (int k = i-1; k <= i+1; ++k) {
                    if(oEnergy[k][j-1].energy == minV){
                        oEnergy[i][j].energy = energy.at<float>(i, j) + minV;

                        oEnergy[i][j].connectedX = k;
                        oEnergy[i][j].connectedY = j-1;
                        break;
                    }
                }
            }
        }
    }

}

bool seam_carving(Mat& in_image, int new_width, int new_height, Mat& out_image){

    // some sanity checks
    // Check 1 -> new_width <= in_image.cols
    if(new_width>in_image.cols){
        cout<<"Invalid request!!! new_width has to be smaller than the current size!"<<endl;
        return false;
    }
    if(new_height>in_image.rows){
        cout<<"Invalid request!!! ne_height has to be smaller than the current size!"<<endl;
        return false;
    }
    
    if(new_width<=0){
        cout<<"Invalid request!!! new_width has to be positive!"<<endl;
        return false;

    }
    
    if(new_height<=0){
        cout<<"Invalid request!!! new_height has to be positive!"<<endl;
        return false;
        
    }

    
    return seam_carving_trivial(in_image, new_width, new_height, out_image);
}


// seam carves by removing trivial seams
bool seam_carving_trivial(Mat& in_image, int new_width, int new_height, Mat& out_image){

    Mat iimage = in_image.clone();
    Mat oimage = in_image.clone();
    Mat gray = Mat(iimage.rows, iimage.cols, CV_32FC1);

    for (int i = 0; i < iimage.rows; ++i) {
        for (int j = 0; j < iimage.cols; ++j) {
            Vec3b pixel = in_image.at<Vec3b>(i, j);
            gray.at<float>(i, j) = cGray(pixel[0], pixel[1], pixel[2]);
        }
    }

    while(iimage.rows!=new_height || iimage.cols!=new_width){
        // horizontal seam if needed
        if(iimage.rows>new_height){
            reduce_horizontal_seam_trivial(gray, iimage, oimage);
            iimage = oimage.clone();
        }
        
        if(iimage.cols>new_width){
            reduce_vertical_seam_trivial(gray, iimage, oimage);
            iimage = oimage.clone();

        }
    }
    
    out_image = oimage.clone();
    return true;
}

// horizontl trivial seam is a seam through the center of the image
bool reduce_horizontal_seam_trivial(Mat &gray, Mat& in_image, Mat& out_image){

    // retrieve the dimensions of the new image
    int rows = in_image.rows;
    int cols = in_image.cols;
    


    Mat energy = Mat(rows, cols, CV_32FC1);
    calculateEnergry(gray, energy, rows, cols);

    Vertex empty = Vertex(-1,-1,-1);
    vector<vector<Vertex > > op_energy(rows, vector<Vertex>(cols, empty));


    optimalHEnergy(energy, op_energy, rows, cols);
    //populate the image

    Vertex sPoint = Vertex(-1,-1,99999);
    for (int i = 0; i < rows; ++i) {
        if(op_energy[i][cols-1].energy < sPoint.energy){
            sPoint = op_energy[i][cols-1];
        }
    }

//    for (int i = 0; i < rows; ++i) {
//        if(i < sPoint.x){
//            out_image.at<Vec3b>(i, cols-1) = in_image.at<Vec3b>(i, cols-1);
//        }else if(i > sPoint.x){
//            out_image.at<Vec3b>(i-1, cols-1) = in_image.at<Vec3b>(i, cols-1);
//            gray.at<float>(i-1, cols-1) = gray.at<float>(i, cols-1);
//        }
//    }
    // create an image slighly smaller
    out_image = Mat(rows-1, cols, CV_8UC3);
    int j = cols-1;
    while(j>=0){

        for (int i = 0; i < rows; ++i) {
            if(i < sPoint.x){
                out_image.at<Vec3b>(i, j) = in_image.at<Vec3b>(i, j);
            }else if(i > sPoint.x){
                out_image.at<Vec3b>(i-1, j) = in_image.at<Vec3b>(i, j);
                gray.at<float>(i-1, j) = gray.at<float>(i, j);
            }
        }
        if(j!=0){
            sPoint = op_energy[sPoint.connectedX][sPoint.connectedY];//sPoint.connected;
        }
        j--;
    }

    return true;
}

// vertical trivial seam is a seam through the center of the image
bool reduce_vertical_seam_trivial(Mat &gray, Mat& in_image, Mat& out_image){
    // retrieve the dimensions of the new image
    int rows = in_image.rows;
    int cols = in_image.cols;

    Mat energy = Mat(rows, cols, CV_32FC1);
    calculateEnergry(gray, energy, rows, cols);

    Vertex empty = Vertex(-1,-1,-1);
    vector<vector<Vertex > > op_energy(rows, vector<Vertex>(cols, empty));


    optimalVEnergy(energy, op_energy, rows, cols);
    Vertex sPoint = Vertex(-1,-1,99999);
    for (int j = 0; j < cols; ++j) {
        if(op_energy[rows-1][j].energy < sPoint.energy){
            sPoint = op_energy[rows-1][j];
        }
    }
    // create an image slighly smaller
    out_image = Mat(rows, cols-1, CV_8UC3);
    int i = rows-1;
    while(i>=0){
        for (int j = 0; j < cols; ++j) {
            if(j<sPoint.y){
                out_image.at<Vec3b>(i, j) = in_image.at<Vec3b>(i, j);
            }else if(j>sPoint.y){
                out_image.at<Vec3b>(i, j-1) = in_image.at<Vec3b>(i, j);
                gray.at<float>(i, j-1) = gray.at<float>(i, j);
            }
        }
        if(i != 0){
            sPoint = op_energy[sPoint.connectedX][sPoint.connectedY];//sPoint = sPoint.connected;
        }

        i--;
    }
    
    return true;
}
