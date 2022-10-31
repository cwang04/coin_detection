#include <stdio.h>
#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;

using namespace cv;

void part1(string filename, int low, int high, double sigma, int spacing){
    int smallR = 83;
    int bigR = 165;
    Mat image;
    image = imread( filename, 1 );
    Mat copy;
    copy = image;

    if ( !image.data )
    {
        printf("No image data \n");
        return;
    }

    imwrite("./copy.jpg", image);
    
    //gray scale image
    Mat gray;
    gray = imread( filename, IMREAD_GRAYSCALE );
    imwrite("./imageg.jpg", gray); 
    
    if(filename=="coins3.jpg"){
        blur(gray, gray, Size(5,5));
    } else if(filename=="coin1.jpg"){
        medianBlur(gray, gray, 5);
        GaussianBlur(gray, gray, Size(5, 5), sigma);
    } else{
        GaussianBlur(gray, gray, Size(5, 5), sigma);
    }
    imwrite("./blurred.jpg", gray); 
    
    //canny edge detection
    Mat canny;
    Canny(gray, canny, low, high, 3);
    imwrite("./imagef.jpg", canny);   
    
    vector<Vec3f> circles;
    //pic 1 sigma = 1.5
    //HoughCircles(gray, circles, HOUGH_GRADIENT, 1, 35, 100, 40, smallR, bigR);
    //pic 2 sigma = 3
    //HoughCircles(gray, circles, HOUGH_GRADIENT, 1, 60, 120, 40, smallR, bigR); 
    //pic 3 
    HoughCircles(gray, circles, HOUGH_GRADIENT, 1, spacing, high, low, smallR, bigR);
    vector<int> coins{0, 0, 0, 0, 0};//silver dollar, quarter, nickel, penny, dime
    double total = 0;
    
    for(int i = 0; i < circles.size(); i++){
        Vec3i c = circles[i];
        Point center = Point(c[0], c[1]);
        int radius = c[2];
        
        if(radius >= 160){//silver dollar
            circle(copy, center, radius, Scalar(0, 255, 255), 3, LINE_AA);
            total+=1.0;
            coins[0]++;
        } else if((radius < 115) && (radius > 103)){//quarter
            circle(copy, center, radius, Scalar(0, 255, 0), 3, LINE_AA);
            total += .25;
            coins[1]++;
        } else if((radius <= 103) && (radius >= 93)){//nickel
            circle(copy, center, radius, Scalar(255, 0, 255), 3, LINE_AA);
            total += .05;
            coins[2]++;
        } else if(radius < 93){
            Mat cropped = image(Range(c[1]-radius, c[1]+radius), Range(c[0]-radius, c[0]+radius));
            Scalar avg;
            avg = mean(cropped);
            if(avg[2]>=(1.1*avg[0])){
                circle(image, center, radius, Scalar(0, 0, 255), 3, LINE_AA);
                total += .01;
                coins[3]++;
            } 
            else{
                circle(image, center, radius, Scalar(255, 0, 0), 3, LINE_AA);
                total += .10;
                coins[4]++;
            }
        }
    }
    cout << coins[3] << " Pennies, " << coins[2] << " Nickels, " << coins[4] << " Dimes, " << coins[1] << " Quarters, " << coins[0];
    cout<<fixed;
    cout<<setprecision(2);
    cout<<" Silver Dollars, Total Sum: $" << total << endl;
    ofstream ofs("output.txt");
    ofs << coins[3] << " Pennies, " << coins[2] << " Nickels, " << coins[4] << " Dimes, " << coins[1] << " Quarters, " << coins[0] << " Silver Dollars, Total Sum: $";
    ofs<<fixed;
    ofs<<setprecision(2);
    ofs << total;
    ofs.close();
    
    imwrite("./coins.jpg", copy);    
    
    waitKey(0);
}

int main(int argc, char** argv )
{  
    if (argc == 2) {
        //lower threshold, upperthreshold, filename, sigma, radius spacing
        part1(argv[1], 40, 100, 1.5, 35);
    }
    else if (argc == 10) {
        int high, low, radius;
        double sigma;
        istringstream iss(argv[3]);
        if (iss >> low) {
        }
        istringstream iss1(argv[5]);
        if (iss1 >> high) {
        }
        istringstream iss2(argv[7]);
        if (iss2 >> radius) {
        }
        istringstream iss3(argv[9]);
        if (iss3 >> sigma) {
        }
        //first 2 values are the upper and lower cannyedge thresholds
        //then it is the filename
        //next is the threshold for the votes for something to be a center
        //then there are the high and low radius integers
        //finally there is the center threshold, which is a percent coverage meaning that it is a double
        part1(argv[1], low, high, sigma, radius);
    } else {
        printf("usage: DisplayImage.out <Image_Path>\n");
        return -1;
    }
    return 0;
}
