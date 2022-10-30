//Colin Wang 
//3/13/22
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <sstream>
#include <chrono>

using namespace std;

class Pixel {
private:
    int r, g, b;
public:
    Pixel(int red, int green, int blue) {
        r = red;
        g = green;
        b = blue;
    }
    Pixel() {
        r = 0;
        g = 0;
        b = 0;
    }
    int getR() {
        return r;
    }
    int getG() {
        return g;
    }
    int getB() {
        return b;
    }
    void setR(int red) {
        r = red;
    }
    void setG(int green) {
        g = green;
    }
    void setB(int blue) {
        b = blue;
    }
};

class Coin {
private:
    int x, y, r, value;
public:
    Coin(int xp, int yp, int rad) {
        x = xp;
        y = yp;
        r = rad;
        value = 0;
    }
    Coin() {
        x = 0;
        y = 0;
        r = 0;
        value = 0;
    }
    int getX() {
        return x;
    }
    int getY() {
        return y;
    }
    int getRadius() {
        return r;
    }
    void setValue(int coin) {
        value = coin;
    }
    int getValue() {
        return value;
    }
};

int rounddouble(double x) {
    return floor(x + 0.5);
}

int sobel_x(vector<int> grid) {//returns the sobel x value
    return (grid[0]) - (grid[2]) + (2 * grid[3]) - (2 * grid[5]) + (grid[6]) - (grid[8]);
}

int sobel_y(vector<int> grid) {//returns the sobel y value
    return (grid[0]) + (2 * grid[1]) + (grid[2]) - (grid[6]) - (2 * grid[7]) - (grid[8]);
}

int grayscale(int r, int g, int b) {//returns grayscale value of rgb
    return (int)((r + b + g) / 3);
}

vector<int> get_grid(int x, int y, int** ppm) {//gets a 3x3 around the x,y point the ppm here already has a 0-border
    vector<int> grid;
    for (int i = y - 1; i < y + 2; i++) {
        for (int j = x - 1; j < x + 2; j++) {
            grid.push_back(ppm[j][i]);
        }
    }
    return grid;
}

int edgeCheck(int threshold, int input) {//returns 1 if its an edge
    if (input > threshold) {
        return 1;
    } return 0;
}

int edgeDetect(int low, int high, int input) {
    if (input > low) {
        if (input < high) {
            return 1;
        } return 2;
    } return 0;
}

void make_ppm(int** arr, string s, int height, int width, int val) {//creates a ppm file given the matrix
    ofstream ofs(s);
    ofs << "P3" << endl << width << ' ' << height << endl << val << endl;

    for (int x = 0; x < height; x++) {
        for (int y = 0; y < width; y++) {
            ofs << arr[x][y] << ' ' << arr[x][y] << ' ' << arr[x][y] << ' ';
        }
        ofs << endl;
    }
    ofs.close();
}

void make_ppm_vector(vector<vector<Pixel>> ppm, string s, int val) {//creates a ppm file given the matrix
    ofstream ofs(s);
    ofs << "P3" << endl << ppm[0].size() << ' ' << ppm.size() << endl << val << endl;

    for (int x = 0; x < ppm.size(); x++) {
        for (int y = 0; y < ppm[0].size(); y++) {
            ofs << ppm[x][y].getR() << ' ' << ppm[x][y].getG() << ' ' << ppm[x][y].getB() << ' ';
        }
        ofs << endl;
    }
    ofs.close();
}

vector<vector<Pixel>> read_file(string filename) {
    ifstream myfile;
    myfile.open(filename);
    string temp;
    myfile >> temp;
    int x, y, r, g, b;
    myfile >> x;
    myfile >> y;
    myfile >> temp;
    vector<vector<Pixel>> ppm(y);
    for (int i = 0; i < y; i++) {
        ppm[i] = vector<Pixel>(x);
        for (int j = 0; j < x; j++) {
            myfile >> r;
            myfile >> g;
            myfile >> b;
            Pixel temp(r, g, b);
            ppm[i][j] = temp;
        }
    }

    return ppm;
}

int** read_gray_file(string filename) {//reads ppm file and returns 2D array of grayscale values
    ifstream myfile;
    myfile.open(filename);
    string temp;
    myfile >> temp;//P3
    int x, y, r, g, b;
    myfile >> x;//width
    myfile >> y;//height also the number of lines
    myfile >> temp;
    int** ppm = new int* [y];
    for (int i = 0; i < y; i++) {
        ppm[i] = new int[x];
        for (int j = 0; j < x; j++) {
            myfile >> r;
            myfile >> g;
            myfile >> b;
            ppm[i][j] = grayscale(r, g, b);
        }
    }

    return ppm;
}

vector<int> get_specifications(string filename) {
    ifstream myfile;
    myfile.open(filename);
    string temp;
    int x, y, val;
    vector<int> dimensions;
    myfile >> temp;
    myfile >> x;
    myfile >> y;
    myfile >> val;
    dimensions.push_back(x);
    dimensions.push_back(y);
    dimensions.push_back(val);
    return dimensions;
}

//Bresenhams Line for Voting
void set_pixel(int i, int j, int** ppm, int height, int width) {
    if ((i >= 0) and (i < height) and (j >= 0) and (j < width)) {
        ppm[i][j] = ppm[i][j] + 1;
    }
}

void draw_line_base(int x1, int x2, int y1, int y2, int** ppm, int height, int width) {//base case of Bresenham's, positive deltaY and deltaX driving
    int j = y1, e = (y2 - y1) - (x2 - x1);
    for (int i = x1; i < x2; i++) {
        set_pixel(i, j, ppm, height, width); //illuminate or fill in the matrix
        if (e >= 0) {
            j += 1;
            e -= (x2 - x1);
        }
        e += (y2 - y1);
    }
}

void draw_line_neg(int x1, int x2, int y1, int y2, int** ppm, int height, int width) { //negative deltaY, deltaX driving axis
    int j = y1, e = (y2 - y1) - (x2 - x1);
    for (int i = x1; i < x2; i++) {
        set_pixel(i, j, ppm, height, width);
        if (e >= 0) {
            j -= 1;
            e -= (x2 - x1);
        }
        e -= (y2 - y1);
    }
}

void draw_line_reverse(int x1, int x2, int y1, int y2, int** ppm, int height, int width) { //positive deltaY, deltaY driving axis
    int j = y1, e = (y2 - y1) - (x2 - x1);
    for (int i = x1; i < x2; i++) {
        set_pixel(j, i, ppm, height, width); //illuminate or fill in the matrix
        if (e >= 0) {
            j += 1;
            e -= (x2 - x1);
        }
        e += (y2 - y1);
    }
}

void draw_line_vert(int x, int y1, int y2, int** ppm, int height, int width) { //vertical case
    for (int i = y1; i < y2; i++) {
        set_pixel(x, i, ppm, height, width);
    }
}

void draw_line_neg_reverse(int x1, int x2, int y1, int y2, int** ppm, int height, int width) {//negative deltaY, deltaY driving axis
    int j = y2, e = (y2 - y1) + (x2 - x1);
    for (int i = x2; i < x1; i++) {
        set_pixel(j, i, ppm, height, width); //illuminate or fill in the matrix
        if (e >= 0) {
            j -= 1;
            e += (x2 - x1);
        }
        e += (y2 - y1);
    }
}

void draw_line(double* p1, double* p2, int** ppm, int height, int width) { //this is the code that we traced over in class
    int x1, x2, y1, y2;
    if (p1[0] < p2[0]) {
        x1 = rounddouble(p1[0] * height), x2 = rounddouble(p2[0] * height), y1 = rounddouble(p1[1] * width), y2 = rounddouble(p2[1] * width);
    }
    else {
        x1 = rounddouble(p2[0] * height), x2 = rounddouble(p1[0] * height), y1 = rounddouble(p2[1] * width), y2 = rounddouble(p1[1] * width);
    }
    set_pixel(x1, y1, ppm, height, width);
    set_pixel(x2, y2, ppm, height, width);
    if (x2 == x1) { //vertical case
        if (y2 > y1) {
            draw_line_vert(x1, y1, y2, ppm, height, width);
        }
        else {
            draw_line_vert(x1, y2, y1, ppm, height, width);
        }
    }
    else if ((y2 - y1) >= 0) {
        if ((x2 - x1) >= (y2 - y1)) {
            draw_line_base(x1, x2, y1, y2, ppm, height, width);
        }
        else {
            draw_line_reverse(y1, y2, x1, x2, ppm, height, width);
        }
    }
    else {
        if (abs(y2 - y1) >= (x2 - x1)) {
            draw_line_neg_reverse(y1, y2, x1, x2, ppm, height, width);
        }
        else {
            draw_line_neg(x1, x2, y1, y2, ppm, height, width);
        }
    }
}

int magnitude(int x, int y) {
    return sqrt(pow(x, 2) + pow(y, 2));
}

int** getMagnitudeGrid(int** border, int height, int width) {//takes the bordered input ppm and returns the magnitude matrix
    int** ppm = new int* [height];
    for (int i = 0; i < height; i++) {
        ppm[i] = new int[width];
        for (int j = 0; j < width; j++) {
            if (i == 0 || i == height - 1 || j == 0 || j == width - 1) {
                ppm[i][j] = 0;
            }
            else {
                ppm[i][j] = magnitude(sobel_x(get_grid(i, j, border)), sobel_y(get_grid(i, j, border)));
            }
        }
    }
    return ppm;
}

int** hysteresis_helper(int** visited, int** ppm, int x, int y, int height, int width) {//y is x and x is y, so y<height and x<width
    if (x < 0 || y < 0 || x >= height || y >= width) {
        return ppm;
    }
    if (visited[x][y] == 1) {
        return ppm;
    }
    if (ppm[x][y] == 0) {
        return ppm;
    }
    else {
        ppm[x][y] = 2;
    }
    visited[x][y] = 1;
    ppm = hysteresis_helper(visited, ppm, x + 1, y, height, width);
    ppm = hysteresis_helper(visited, ppm, x + 1, y + 1, height, width);
    ppm = hysteresis_helper(visited, ppm, x + 1, y - 1, height, width);
    ppm = hysteresis_helper(visited, ppm, x - 1, y, height, width);
    ppm = hysteresis_helper(visited, ppm, x - 1, y + 1, height, width);
    ppm = hysteresis_helper(visited, ppm, x - 1, y - 1, height, width);
    ppm = hysteresis_helper(visited, ppm, x, y + 1, height, width);
    ppm = hysteresis_helper(visited, ppm, x, y - 1, height, width);
    return ppm;
}

int** hysteresis(int** ppm, int height, int width) {
    int** output = new int* [height];
    int** visited = new int* [height];
    for (int x = 0; x < height; x++) {
        visited[x] = new int[width];
        output[x] = new int[width];
        for (int y = 0; y < width; y++) {
            output[x][y] = ppm[x][y];
            visited[x][y] = 0;
        }
    }
    for (int x = 0; x < height; x++) {
        for (int y = 0; y < width; y++) {
            if (visited[x][y] == 0 && output[x][y] == 2) {
                output = hysteresis_helper(visited, output, x, y, height, width);
            }
        }
    }
    for(int i = 0; i < height; i++){
        delete[] visited[i];
    }
    delete[] visited;
    return output;
}

int promote(int input) {//makes the 0, 1, 2 into just 0 and 1 with 1s being converted to 
    if (input > 0) {
        return input - 1;
    } return 0;
}

int roundDegree(double angle) {
    return round(angle / 45) * 45;
}

double** getDegreeGrid(int** ppm, int height, int width) {
    double** output = new double* [height];
    for (int i = 0; i < height; i++) {
        output[i] = new double[width];
        for (int j = 0; j < width; j++) {
            if (i == 0 || j == 0 || i == height - 1 || j == width - 1) {
                output[i][j] = 0;
            }
            else {
                output[i][j] = atan2(sobel_y(get_grid(i, j, ppm)), sobel_x(get_grid(i, j, ppm))) * 180.0 / M_PI;
            }
        }
    }
    return output;
}

int** roundDegreeGrid(double** degree, int height, int width) {
    int** output = new int* [height];
    for (int i = 0; i < height; i++) {
        output[i] = new int[width];
        for (int j = 0; j < width; j++) {
            output[i][j] = roundDegree(degree[i][j]);
        }
    }
    return output;
}

int** nonmaxSupression(int** degree, int** magnitude, int height, int width) {
    int** nms = new int* [height];
    for (int i = 0; i < height; i++) {
        nms[i] = new int[width];
        for (int j = 0; j < width; j++) {
            if (i == 0 || j == 0 || i == height - 1 || j == width - 1) {
                nms[i][j] = 0;
            }
            else if (degree[i][j] == 0 || degree[i][j] == 180 || degree[i][j] == -180) { //0 and 180 should be the one left and right [i][j-1] and [i][j+1]
                if ((magnitude[i][j] >= magnitude[i][j - 1]) && (magnitude[i][j] >= magnitude[i][j + 1])) {
                    nms[i][j] = magnitude[i][j];
                }
                else {
                    nms[i][j] = 0;
                }
            }
            else if (degree[i][j] == 45 || degree[i][j] == -135) {//45 and -135 should be the one diagonally right so [i-1][j-1] and [i+1][j+1]
                if ((magnitude[i][j] >= magnitude[i - 1][j - 1]) && (magnitude[i][j] >= magnitude[i + 1][j + 1])) {
                    nms[i][j] = magnitude[i][j];
                }
                else {
                    nms[i][j] = 0;
                }
            }
            else if (degree[i][j] == -45 || degree[i][j] == 135) {//-45 and 135 should be diagonally left so [i-1][j+1] and [i+1][j-1]
                if ((magnitude[i][j] >= magnitude[i - 1][j + 1]) && (magnitude[i][j] >= magnitude[i + 1][j - 1])) {
                    nms[i][j] = magnitude[i][j];
                }
                else {
                    nms[i][j] = 0;
                }
            }
            else {// 90 and -90 so [i+1][j] and [i-1][j]
                if ((magnitude[i][j] >= magnitude[i - 1][j]) && (magnitude[i][j] >= magnitude[i + 1][j])) {
                    nms[i][j] = magnitude[i][j];
                }
                else {
                    nms[i][j] = 0;
                }
            }
        }
    }
    return nms;
}

int** combine(int** nms, int** hysteresis, int height, int width) {
    int** finish = new int* [height];
    for (int i = 0; i < height; i++) {
        finish[i] = new int[width];
        for (int j = 0; j < width; j++) {
            finish[i][j] = (nms[i][j] + hysteresis[i][j]) / 2;
        }
    }
    return finish;
}

//x and y coordinates of the cell in the matrix, in double form for calculations
void getLine(int x, int y, double** degrees, int** voting, int height, int width) {
    double slope = tan(degrees[x][y] * M_PI / 180.0);
    //cout<<slope<<endl;
    double c = ((double)y / (double)width) - (slope * ((double)x / (double)height));
    double ymin = c;
    double ymax = c + slope;
    double point1[2] = { 0.0, ymin };
    double point2[2] = { 1.0, ymax };
    draw_line(point1, point2, voting, height, width);
}

int** votingGrid(double** degree, int height, int width, int** ppm) {
    int** voting = new int* [height];
    for (int i = 0; i < height; i++) {
        voting[i] = new int[width];
        for (int j = 0; j < width; j++) {
            voting[i][j] = 0;
        }
    }
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            if (ppm[i][j] == 1) {
                getLine(i, j, degree, voting, height, width);
            }
        }
    }
    return voting;
}

void draw_pixel(int i, int j, vector<vector<Pixel>>& ppm, int color[]) {
    if ((i >= 0) and (i < ppm.size()) and (j >= 0) and (j < ppm[0].size())) {
        Pixel temp(color[0], color[1], color[2]);
        ppm[i][j] = temp;
    }
}

void draw_circle(int xp, int yp, int radius, vector<vector<Pixel>>& ppm, int color[]) {
    int xmax = (int)radius / sqrt(2);
    int y = radius;
    int y2 = y * y;
    int ty = (2 * y) - 1;
    int y2_new = y2;

    for (int x = 0; x <= xmax + 2; x++) {
        if ((y2 - y2_new) >= ty) {
            y2 -= ty;
            y -= 1;
            ty -= 2;
        }
        draw_pixel(xp + x, yp + y, ppm, color);
        draw_pixel(xp + x, yp - y, ppm, color);
        draw_pixel(xp - x, yp + y, ppm, color);
        draw_pixel(xp - x, yp - y, ppm, color);
        draw_pixel(xp + y, yp + x, ppm, color);
        draw_pixel(xp + y, yp - x, ppm, color);
        draw_pixel(xp - y, yp + x, ppm, color);
        draw_pixel(xp - y, yp - x, ppm, color);

        y2_new -= (2 * x) - 3;
    }
}

bool checkCenter(int x, int y, int** ppm, int height, int width) {
    if (x >= 0 && y >= 0 && x < height && y < width) {
        if (ppm[x][y] == 1) {
            return true;
        }
    } return false;
}

double circle_count(int xp, int yp, int radius, int** ppm, int height, int width) {
    int xmax = (int)radius / sqrt(2);
    int y = radius;
    int y2 = y * y;
    int ty = (2 * y) - 1;
    int y2_new = y2;
    int count = 0;
    int total = 0;
    for (int x = 0; x <= xmax + 2; x++) {
        if ((y2 - y2_new) >= ty) {
            y2 -= ty;
            y -= 1;
            ty -= 2;
        }
        if (checkCenter(xp + x, yp + y, ppm, height, width)) {
            count++;
        }
        if (checkCenter(xp + x, yp - y, ppm, height, width)) {
            count++;
        }
        if (checkCenter(xp - x, yp + y, ppm, height, width)) {
            count++;
        }
        if (checkCenter(xp - x, yp - y, ppm, height, width)) {
            count++;
        }
        if (checkCenter(xp + y, yp + x, ppm, height, width)) {
            count++;
        }
        if (checkCenter(xp + y, yp - x, ppm, height, width)) {
            count++;
        }
        if (checkCenter(xp - y, yp + x, ppm, height, width)) {
            count++;
        }
        if (checkCenter(xp - y, yp - x, ppm, height, width)) {
            count++;
        }
        total += 8;
        y2_new -= (2 * x) - 3;
    }
    //cout<<radius<<": "<<total<<", "<<(double)count/(double)total<<endl;
    return (double)count / (double)total;
}

int** findCenters(int** voting, int height, int width, int threshold) {
    int** center = new int* [height];
    for (int i = 0; i < height; i++) {
        center[i] = new int[width];
        for (int j = 0; j < width; j++) {
            if (voting[i][j] > threshold) {
                center[i][j] = 1;
            }
            else {
                center[i][j] = 0;
            }
        }
    }
    return center;
}

int** blur(int** ppm, int height, int width) {
    double** kernel = new double* [5];
    for (int x = 0; x < 5; x++) {
        kernel[x] = new double[5];
    }
    int** output = new int* [height];
    for (int x = 0; x < height; x++) {
        output[x] = new int[width];
        for (int y = 0; y < width; y++) {
            output[x][y] = ppm[x][y];
        }
    }
    double sigma = 1.5;
    double sum = 0.0, q = 2.0 * sigma * sigma;
    double p;
    for (int x = -2; x < 3; x++) {
        for (int y = -2; y < 3; y++) {
            p = sqrt(x * x + y * y);
            kernel[2 + x][2 + y] = (exp(-(p * p) / q)) / (M_PI * q);
            sum += kernel[2 + x][2 + y];
        }
    }
    for (int x = 0; x < 5; x++) {
        for (int y = 0; y < 5; y++) {
            kernel[x][y] /= sum;
        }
    }
    for (int i = 2; i < height - 2; i++) {
        for (int j = 2; j < width - 2; j++) {
            sum = 0.0;
            for (int x = -2; x < 3; x++) {
                for (int y = -2; y < 3; y++) {
                    sum += (kernel[2 + x][2 + y] * ppm[i + x][y + j]);
                }
            }
            output[i][j] = (int)sum;
        }
    }
    for(int i = 0; i < 5; i++){
        delete[] kernel[i];
    }
    delete[] kernel;
    return output;
}

vector<Coin> remove_duplicates(vector<Coin> ppm, double*** voting) {
    vector<Coin> output;
    int xmax, ymax, tempx, tempy, prevxmax, prevymax;
    for (int x = 0; x < ppm.size(); x++) {
        if (output.size() == -0) {
            output.push_back(Coin(ppm[x].getX(), ppm[x].getY(), ppm[x].getRadius()));
        }
        else if ((abs(output[output.size() - 1].getX() - ppm[x].getX()) > 40) && (abs(output[output.size() - 1].getY() - ppm[x].getY())) > 40) {
            output.push_back(Coin(ppm[x].getX(), ppm[x].getY(), ppm[x].getRadius()));
        }
    }
    return output;
}

void part2(int high, int low, string filename, int threshold, int bigR, int smallR, double circlethreshold) {//detect circle radius as well as coins based on pixel color and size
    vector<int> specifications = get_specifications(filename);
    int red[3] = { 255, 0, 0 };
    int yellow[3] = { 255, 255, 0 };
    int purple[3] = { 255, 0, 255 };
    int green[3] = { 0, 255, 0 };
    int blue[3] = { 0, 0, 255 };
    int height = specifications[1];
    int width = specifications[0];//nested array length
    int rgb = specifications[2];
    int** ppm = read_gray_file(filename);//creates matrix of grayscale values of 
    make_ppm(ppm, "imageg.ppm", height, width, rgb);
    int** blurred = blur(ppm, height, width);
    blurred = blur(blurred, height, width);
    make_ppm(blurred, "blur.ppm", height, width, rgb);
    int** mag = getMagnitudeGrid(blurred, height, width);//creates a matrix with the magnitude values
    double** fullDegree = getDegreeGrid(blurred, height, width);//creates a matrix with the degree values
    int** degree = roundDegreeGrid(fullDegree, height, width);//rounds these values to the closest multiple of 45
    //nms grid
    int** nms = nonmaxSupression(degree, mag, height, width);
    make_ppm(nms, "image2.ppm", height, width, 1);

    //generate the hysteresis grid
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            mag[i][j] = edgeDetect(low, high, nms[i][j]);
        }
    }
    int** hysterized = hysteresis(mag, height, width);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            hysterized[i][j] = promote(hysterized[i][j]);
        }
    }
    make_ppm(hysterized, "image1.ppm", height, width, 1);

    //combine the hysteresis and the nms grids
    make_ppm(hysterized, "imagef.ppm", height, width, 1);
    int** voting = votingGrid(fullDegree, height, width, hysterized);
    int max = -1;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (voting[i][j] > max) {
                max = voting[i][j];
            }
        }
    }
    make_ppm(voting, "imagev.ppm", height, width, max);
    int** centers = findCenters(voting, height, width, threshold);

    vector<vector<Pixel>> original = read_file(filename);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (centers[i][j] == 1) {
                for (int x = 1; x <= 4; x++) {
                    draw_circle(i, j, x, original, red);
                }
            }
        }
    }
    make_ppm_vector(original, "imageCC.ppm", rgb);

    //part2
    double*** voting2 = new double** [height];
    vector<Coin> centers2;
    double tempmax;
    for (int i = 0; i < height; i++) {
        voting2[i] = new double* [width];
        for (int j = 0; j < width; j++) {
            voting2[i][j] = new double[bigR - smallR + 1];
            if (centers[i][j] == 1) {
                tempmax = -1.0;
                for (int k = 0; k < bigR - smallR; k++) {
                    voting2[i][j][k] = circle_count(i, j, smallR + k, combined, height, width);//returns percent of circle that is an edge
                    if (voting2[i][j][k] > tempmax) {
                        tempmax = voting2[i][j][k];
                        if (tempmax > circlethreshold) {
                            centers2.push_back(Coin(i, j, k));
                        }
                    }
                }
            }
        }
    }
    vector<Coin> cleaned = remove_duplicates(centers2, voting2);
    original = read_file(filename);
    //remove duplicate centers

    //pennies: r~=85 / nickels: r~=95 / quarters: r ~= 108 / dimes: r ~= 80 / half-dollar: r~=136
    for (int x = 0; x < cleaned.size(); x++) {
        int radius = cleaned[x].getRadius() + smallR;
        if (radius >= 125) {
            cleaned[x].setValue(100);//halfdollar
        }
        else if (radius < 125 && radius > 104) {
            cleaned[x].setValue(25);//quarter
        }
        else if (radius <= 104 && radius > 94) {
            cleaned[x].setValue(5);//nickel
        }
        else {
            for (int i = 0; i < 5; i++) {
                if ((original[cleaned[x].getX() + i][cleaned[x].getY() + i].getR() - original[cleaned[x].getX() + i][cleaned[x].getY() + i].getB()) >= 17) {//pennies tend to have large R small B
                    cleaned[x].setValue(1);//penny
                }
                else {
                    cleaned[x].setValue(10);//dime
                }
            }
        }
    }

    //draw coin outlines red for penny, purple for nickel, blue for dime, green for quarters and yellow for silver dollar.
    //additionally create a count for each coin and a running total
    int coins[5] = { 0, 0, 0, 0, 0 };
    double total = 0.0;
    for (int x = 0; x < cleaned.size(); x++) {
        if (cleaned[x].getValue() == 1) {//penny
            draw_circle(cleaned[x].getX(), cleaned[x].getY(), cleaned[x].getRadius() + smallR, original, red);
            total += .01;
            coins[0]++;
        }
        else if (cleaned[x].getValue() == 5) {//nickel
            draw_circle(cleaned[x].getX(), cleaned[x].getY(), cleaned[x].getRadius() + smallR, original, purple);
            total += .05;
            coins[1]++;
        }
        else if (cleaned[x].getValue() == 1) {//dime
            draw_circle(cleaned[x].getX(), cleaned[x].getY(), cleaned[x].getRadius() + smallR, original, blue);
            total += .10;
            coins[2]++;
        }
        else if (cleaned[x].getValue() == 25) {//quarter
            draw_circle(cleaned[x].getX(), cleaned[x].getY(), cleaned[x].getRadius() + smallR, original, green);
            total += .25;
            coins[3]++;
        }
        else if (cleaned[x].getValue() == 100) {//silver dollar
            draw_circle(cleaned[x].getX(), cleaned[x].getY(), cleaned[x].getRadius() + smallR, original, yellow);
            total += 1.0;
            coins[4]++;
        }
    }
    make_ppm_vector(original, "coins.ppm", rgb);
    ofstream ofs("output.txt");
    ofs << coins[0] << " Pennies, " << coins[1] << " Nickels, " << coins[2] << " Dimes, " << coins[3] << " Quarters, " << coins[4] << " Silver Dollars, Total Sum: $" << total;
    ofs.close();
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        //(int high, int low, string filename, int threshold, int bigR, int smallR, double circlethreshold)
        //pic1: 100, 50, "coin1.ppm", 60, 140, 75, .25
        //pic2: 
        part2(130, 32, "coin1.ppm", 66, 140, 85, .15);
    }
    else if (argc == 15) {
        int high, low, threshold, highR, lowR;
        double center;
        istringstream iss(argv[2]);
        if (iss >> low) {
        }
        istringstream iss1(argv[4]);
        if (iss1 >> high) {
        }
        istringstream iss2(argv[8]);
        if (iss2 >> threshold) {
        }
        istringstream iss3(argv[10]);
        if (iss2 >> highR) {
        }
        istringstream iss4(argv[12]);
        if (iss2 >> lowR) {
        }
        istringstream iss5(argv[14]);
        if (iss2 >> center) {
        }
        //first 2 values are the upper and lower cannyedge thresholds
        //then it is the filename
        //next is the threshold for the votes for something to be a center
        //then there are the high and low radius integers
        //finally there is the center threshold, which is a percent coverage meaning that it is a double
        part2(high, low, argv[6], threshold, highR, lowR, center);
    }
}
