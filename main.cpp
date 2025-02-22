#include <iostream>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <chrono>
#include <string>
#include <thread>
#include <queue>

using vector2D = std::vector< std::vector<int> >;
using string = std::string;

struct matrixSandpile{
    std::vector<int> verticies;
    std::vector< std::vector<int> > congruence_map;
    std::vector<int> odometer;
};

// Overloads

std::vector<int> operator+(std::vector<int> a, std::vector<int> b){
    for (int i = 0; i < a.size(); i++) {
        a[i] += b[i];
    };
    return a;
};

std::vector<int> operator-(std::vector<int> a, std::vector<int> b){
    for (int i = 0; i < a.size(); i++) {
        a[i] -= b[i];
    };
    return a;
};

// Helper Functions

void createBMP(const vector2D imageData, const vector2D colours, const string filename) {
    int height = imageData.size();
    int width = imageData[0].size();

    // BMP file header
    uint8_t bmpFileHeader[14] = {
        0x42, 0x4D,             // Signature "BM"
        0, 0, 0, 0,             // File size in bytes (placeholder)
        0, 0,                   // Reserved
        0, 0,                   // Reserved
        54, 0, 0, 0             // Offset to pixel data (54 bytes)
    };

    // BMP info header (40 bytes)
    uint8_t bmpInfoHeader[40] = {
        40, 0, 0, 0,            // Info header size
        0, 0, 0, 0,             // Image width
        0, 0, 0, 0,             // Image height
        1, 0,                   // Planes
        24, 0,                  // Bits per pixel
        0, 0, 0, 0,             // Compression
        0, 0, 0, 0,             // Image size
        0x13, 0x0B, 0, 0,       // Horizontal resolution
        0x13, 0x0B, 0, 0,       // Vertical resolution
        0, 0, 0, 0,             // Colors in palette
        0, 0, 0, 0              // Important colors
    };

    // Write correct width and height to header
    bmpInfoHeader[4] = width & 0xFF;
    bmpInfoHeader[5] = (width >> 8) & 0xFF;
    bmpInfoHeader[6] = (width >> 16) & 0xFF;
    bmpInfoHeader[7] = (width >> 24) & 0xFF;

    bmpInfoHeader[8] = height & 0xFF;
    bmpInfoHeader[9] = (height >> 8) & 0xFF;
    bmpInfoHeader[10] = (height >> 16) & 0xFF;
    bmpInfoHeader[11] = (height >> 24) & 0xFF;


    int rowPadding = (4 - (width * 3) % 4) % 4;
    int fileSize = 54 + (width * 3 + rowPadding) * height;

    // Write correct file size to header
    bmpFileHeader[2] = fileSize & 0xFF;
    bmpFileHeader[3] = (fileSize >> 8) & 0xFF;
    bmpFileHeader[4] = (fileSize >> 16) & 0xFF;
    bmpFileHeader[5] = (fileSize >> 24) & 0xFF;

    // Open file and write headers
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        std::cerr << "Error: Cannot open file for writing: " << filename << std::endl;
        return;
    }
    outFile.write(reinterpret_cast<char*>(bmpFileHeader), sizeof(bmpFileHeader));
    outFile.write(reinterpret_cast<char*>(bmpInfoHeader), sizeof(bmpInfoHeader));

    // Write mapped pixel data
    // BMP goes bottom up
    for (int i = height - 1; i >= 0; --i) {
        for (int j = 0; j < width; ++j) {
            // RGB individually
            outFile.put(colours[imageData[i][j]][0]);
            outFile.put(colours[imageData[i][j]][1]);
            outFile.put(colours[imageData[i][j]][2]);
        }
        // Pad
        for (int k = 0; k < rowPadding; ++k) {
            outFile.put(0);
        }
    }

    outFile.close();
    std::cout << "BMP file created successfully: " << filename << std::endl;
}

vector2D readBMP(const string &filename, const vector2D colours) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file!" << std::endl;
        exit(1);
    }

    // Read BMP Header (First 54 bytes)
    unsigned char header[54];
    file.read(reinterpret_cast<char *>(header), 54);

    // Extract width and height (little-endian format)
    int width = header[18] | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
    int height = header[22] | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);

    // Each row is padded to a multiple of 4 bytes
    int row_padded = (3 * width + 3) & (~3);

    // Store pixel data
    vector2D pixels;
    std::vector<int> blank_row(width, 0);
    for (int i=0; i < height; i++){
        pixels.push_back(blank_row);
    }
    std::vector<int> pix(3, 0);

    for (int i = 0; i < height; i++) {
        std::vector<unsigned char> row(row_padded);
        file.read(reinterpret_cast<char *>(row.data()), row_padded);

        for (int j = 0; j < width; j++) {
            pix[2] = row[j * 3];
            pix[1] = row[j * 3 + 1];
            pix[0] = row[j * 3 + 2];
            pixels[height - 1 - i][j] = std::find(colours.begin(), colours.end(), pix) - colours.begin(); // Flip vertically
        }
    }

    file.close();
    return pixels;
}

void createTXT(const std::vector<int> data, const string filename){
    std::ofstream outFile(filename, std::ios::binary);
    for (int a : data) {
        outFile << a << ",";
    }
    outFile.close();
}

std::vector<int> readTXT(const string filename){
    std::vector<int> numbers;
    std::ifstream file(filename);
    
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return numbers;
    }

    std::string line, value;
    while (std::getline(file, line)) { // Read the whole line
        std::stringstream ss(line);
        while (std::getline(ss, value, ',')) { // Split by comma
            numbers.push_back(std::stoi(value)); // Convert to int and store
        }
    }

    return numbers;
}

vector2D nColourScheme(const int colours){
    vector2D scheme;
    for (int level = 1; level < colours+1; level++){
        std::vector<int> map;
        map.push_back(255/level);
        map.push_back(255/level);
        map.push_back(255/level);
        scheme.push_back(map);
    };
    return scheme;
};

vector2D unflatten(const std::vector<int> sandpile, const int width, const int height){
    std::vector<int> row = {sandpile[0]};
    vector2D grid;

    for (int i = 1; i < width*height; i++){
        // If the vertex is the start of a new row, add the current row to the gird
        // and start a new row with that vertex as the first element
        // Otherwise add the vertex to the current row
        if (i % width == 0){
            grid.push_back(row);
            row = {sandpile[i]};
        }
        else {
            row.push_back(sandpile[i]);
        };
    }
    grid.push_back(row);

    return grid;
};

std::vector<int> flatten(const vector2D grid){
    std::vector<int> flattened_grid;
    for (int y = 0; y < grid.size(); y++) {
        for (int x = 0; x < grid[0].size(); x++) {
            flattened_grid.push_back(grid[y][x]);
        }
    }
    return flattened_grid;
}

std::vector<int> scaleOdometer(std::vector<int> odometer, int width, int height, int scale) {
    std::vector<int> scaled((width*scale)*(height*scale), 0);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int s = 0; s < scale; s++) {
                for (int ss = 0; ss < scale; ss++) {
                    scaled[(x * scale)+ ((y+ss) * width * scale) + s] = odometer[(y*width) + x];
                }
            }
        }
    }

    return scaled;
}

std::vector<int> increaseOdometer(std::vector<int> odometer, int width, int height, int increase) {
    std::vector<int> scaled((width+increase)*(height+increase), 0);

    for (int a = 0; a < height; a++) {
        for (int b = 0; b < width; b++) {
            scaled[(width+increase)*(a+1) + 1 + b] = odometer[a*(width) + b];
        }
    }

    return scaled;
}

// Congruence map generators

vector2D generateCongruenceMap(vector2D laplacian){
    vector2D congruence;
    std::vector<int> blank;
    for (int i = 0; i < laplacian.size(); i++){
        congruence.push_back(blank);
        for (int ii = 0; ii < laplacian[i].size(); ii++){
            if (laplacian[i][ii] > 0) {
                congruence[i].push_back(ii);
            }
        };
    };
    return congruence;
};

vector2D generateSquareCongruenceMap(int width, int height) {
    vector2D congruence;
    std::vector<int> blank;

    for (int vert = 0; vert < height * width; vert++){
        congruence.push_back(blank);
        int vm = vert % width;
        int vf = vert / width;
        if (vm != width - 1){
            congruence[vert].push_back(vert+1);
        }
        if (vm != 0){
            congruence[vert].push_back(vert-1);
        }
        if (vf != height - 1){
            congruence[vert].push_back(vert+width);
        }
        if (vf != 0){
            congruence[vert].push_back(vert-width);
        }
    };

    return congruence;
};

vector2D generateTorusCongruenceMap(int width, int height, int offset = 0) {
    vector2D congruence;
    std::vector<int> blank;

    for (int vert = 0; vert < width*height; vert++){
        congruence.push_back(blank);
        if (vert == 0) {
            continue;
        }

        int vm = vert % width;
        int vf = vert / width;
        // Add vertex one to the right
        congruence[vert].push_back( ((vm + 1) % width) + (vf * width));
        // Add vertex one to the left
        if (vm == 0){
            congruence[vert].push_back( width - 1 + (vf * width));
        } else {
            congruence[vert].push_back( ((vm - 1) % width) + (vf * width));
        }
        // Add vertex one down
        if (vm == 0) {
            congruence[vert].push_back( (((vf + 1) % height) * width));
        } else {
            congruence[vert].push_back( ((vm - offset) % width) + (((vf + 1) % height) * width));
        }
        // Add vertex one up
        if (vf == 0) {
            congruence[vert].push_back( ((vm + offset) % width) + ((height - 1) * width));
        } else {
            congruence[vert].push_back( ((vm + offset) % width) + (((vf - 1) % height) * width));
        }
    };

    return congruence;
};

vector2D generateDirectedTorusCongruenceMap(int size, int total_size, int offset) {
    vector2D congruence;
    std::vector<int> blank;

    for (int vert = 0; vert < total_size; vert++){
        congruence.push_back(blank);
        if (vert == 0) {
            continue;
        }

        int vm = vert % size;
        int vf = vert / size;
        
        congruence[vert].push_back((vm * size) + ((vf + 1)%size));
        congruence[vert].push_back((((vm+1) % size) * size) + vf);
    };

    return congruence;
};

vector2D generateBinaryHyperbolicMap(int layers){
    /* // Cmax for binary hyperbolic
   double_cmax[0] = 4;

    for (int i = 1; i < size; i++){
        double_cmax[pow(2, i) - 1] = 8;
        double_cmax[pow(2, i+1) - 2] = 8;
    }

    for (int i = pow(2, size-1)-1; i < pow(2, size); i++){
        double_cmax[i] = 6;
    } */
    vector2D congruence;
    std::vector<int> blank;

    for (int i=0; i < layers; i++){
        for (int vertex=pow(2, i); vertex < pow(2, i+1); vertex++){
            int cvertex = vertex - 1;
            congruence.push_back(blank);
            if (vertex != pow(2, i)) {
                congruence[cvertex].push_back(cvertex - 1);
            }
            if (vertex != pow(2, i+1) - 1) {
                congruence[cvertex].push_back(cvertex + 1);
            }
            if (i != 0) {
                congruence[cvertex].push_back(vertex / 2 - 1);
            }
            if (i != layers - 1) {
                congruence[cvertex].push_back(vertex * 2);
                congruence[cvertex].push_back(vertex * 2 - 1);
            }
        };
    };
    return congruence;
};

vector2D generateCylinderCongruenceMap(int size, int total_size) {
     vector2D congruence;
    std::vector<int> blank;

    for (int vert = 0; vert < total_size; vert++){
        congruence.push_back(blank);
        if (vert == 0){
            continue;
        }

        int vm = vert % size;
        int vf = vert / size;
        if (vm != size - 1){
            congruence[vert].push_back(vert+1);
        } else {
            congruence[vert].push_back(vf * size);
        };
        if (vm != 0){
            congruence[vert].push_back(vert-1);
        } else {
            congruence[vert].push_back(vert + size - 1);
        };
        if (vf != size - 1){
            congruence[vert].push_back(vert+size);
        };
        if (vf != 0){
            congruence[vert].push_back(vert-size);
        };
    };

    return congruence;
}

vector2D generateMobiusCongruenceMap(int size, int total_size) {
     vector2D congruence;
    std::vector<int> blank;

    for (int vert = 0; vert < total_size; vert++){
        congruence.push_back(blank);
        if (vert == 0){
            continue;
        }

        int vm = vert % size;
        int vf = vert / size;
        if (vm != size - 1){
            congruence[vert].push_back(vert+1);
        } else {
            congruence[vert].push_back((size-vf-1) * size);
        };
        if (vm != 0){
            congruence[vert].push_back(vert-1);
        } else {
            congruence[vert].push_back((size - vf)*size - 1);
        };
        if (vf != size - 1){
            congruence[vert].push_back(vert+size);
        };
        if (vf != 0){
            congruence[vert].push_back(vert-size);
        };
    };

    return congruence;
}

// Big boy stuff

void toppleWorker(  matrixSandpile &sandpile, const std::vector<int> &search_range, int &tcount, int &runthru,
                    std::vector<std::mutex> &vertical_locks, std::vector<std::mutex> &horizontal_locks,
                    const int width, const int height, const int search_width, const int search_height, const int n){
    int vm, vf, vmm, vfm;
    bool ver_locking, hor_locking;

    for (int vertex : search_range){
        if (sandpile.verticies[vertex] > 3){
            tcount += 1;
            runthru = 1;

            vm = (vertex % width);
            vf = (vertex / width);

            if (vf == 0 || vf == height - 1 || vm == 0 || vm == width - 1) {
                ver_locking = false;
                hor_locking = false;
            } else {
                vmm = vm  % (search_width);
                vfm = vf  % (search_height);
                ver_locking = (vmm == 0 || vmm == search_width - 1 || vmm == 1 || vmm == search_width - 2);
                hor_locking = (vfm == 0 || vfm == search_height - 1 || vfm == 1 || vfm == search_height - 2);
            }

            if (ver_locking) {
                //std::cout << "locking " << (vmm > 2) + vm / search_width << '\n';
                vertical_locks[(vmm > 2) + (vm / search_width)].lock();
            }
            if (hor_locking) {
                //std::cout << "locking " << (vfm > 2) + vf / search_height << '\n';
                horizontal_locks[(vfm > 2) + (vf / search_height)].lock();
            }
            //
            
            sandpile.odometer[vertex] += 1;
            sandpile.verticies[vertex] -= 4;
            // if ((vertex / size == 0) || (vertex / size == size - 1)) {
            //    sandpile.verticies[vertex] -= 3;
            //} else {
            // };

            for (int topples : sandpile.congruence_map[vertex]){
                sandpile.verticies[topples] += 1;
            };

            if (hor_locking) {
                //std::cout << "unlocking " << (vfm > 2) + vf / search_height << '\n';
                horizontal_locks[(vfm > 2) + (vf / search_height)].unlock();
            }

            if (ver_locking) {
                //std::cout << "unlocking " << (vmm > 2) + vm / search_width << '\n';
                vertical_locks[(vmm > 2) + (vm / search_width)].unlock();
            }
        };
    }
};

std::vector<int> generateToppleRange(int starting_point, int check_width, int check_height, int width){
    std::vector<int> range;
    for (int offset = 0; offset < check_height * check_width; offset++) {
        range.push_back(starting_point + (offset % check_width) + (offset / check_width)*width);
    }
    return range;
};

vector2D generateGridCover(const int covers, const int width, const int height) {
    // generates cover^2 non overlapping rectangles
    vector2D cover;
    int search_width;
    int search_height;

    for (int y = 0; y < covers; y++){
        if (y == covers - 1){ search_height = height - (height / covers * (covers - 1));} else { search_height = height / covers;}
        for (int x = 0; x < covers; x++){
            if (x == covers - 1) { search_width = width - (width / covers * (covers - 1));} else { search_width = width / covers;}
            cover.push_back(generateToppleRange((width * y * (height / covers)) + (x * (width / covers)), search_width, search_height, width));
        }
    }
    
    return cover;
};

matrixSandpile stabaliseMatrixSandpile( matrixSandpile sandpile, const int width, const int height,
                                        const bool pre_topple, const bool debug = false){
    
    
    const int COVERS = 1;
    int runthru = 1;
    int lcount = 0;
    int tcount = 0;
    int total_size = width*height;

    if (debug) {printf("Stabalizing sandpile.\n");}
    if (pre_topple){
        if (debug) {printf("Pre toppling vericies....\n");}
        for (int vertex = 0; vertex < total_size; vertex++) {
            sandpile.verticies[vertex] -= 4*sandpile.odometer[vertex];
            for (int topples : sandpile.congruence_map[vertex]){
                sandpile.verticies[topples] += sandpile.odometer[vertex];
            }
        }
        if (debug) {printf("Done!\n");}
    }

    if (debug) {printf("Toppling vertexes...\n");}

    std::vector<std::thread> threads;
    vector2D gridCover = generateGridCover(COVERS, width, height);

    std::vector<std::mutex> vertical_locks(std::max(2, COVERS*COVERS));
    std::vector<std::mutex> horizontal_locks(std::max(2, COVERS*COVERS));

    while (runthru) {
        lcount += 1;
        runthru = 0;

        for (std::vector<int> &q : gridCover) {
            threads.push_back(std::thread(  &toppleWorker, std::ref(sandpile), std::ref(q), std::ref(tcount), std::ref(runthru),
                                            std::ref(vertical_locks), std::ref(horizontal_locks),
                                            width, height, width / COVERS, height / COVERS, COVERS));
        }

        for (auto&& thread : threads) {
            thread.join();
        }

        threads.clear();
    }
    if (debug) {printf("Done! Sandpile stabalizsed. %d loops perforemd with %d topples.\n", lcount, tcount);}
    return sandpile;
};

std::vector<int> generateNeutralElement(const vector2D cong_map, const int width, const int height, const bool estimation,
                                        std::vector<int>& firstOdo, std::vector<int>& secondOdo, const bool debug = false){
    std::vector<int> double_cmax(width*height, 8);
    std::vector<int> zero_odo(width*height,0);
    matrixSandpile firstStab;
    matrixSandpile secondStab;
    int total_size = width * height;

    // Uncomment for a sink at 0
    double_cmax[0] = 0;

    if (estimation){
        firstStab = stabaliseMatrixSandpile({double_cmax, cong_map, firstOdo}, width, height, true, debug);
        secondStab = stabaliseMatrixSandpile({double_cmax - firstStab.verticies, cong_map, secondOdo}, width, height, true, debug);
    } else {
        firstStab = stabaliseMatrixSandpile({double_cmax, cong_map, zero_odo}, width, height, false, debug);
        secondStab = stabaliseMatrixSandpile({double_cmax - firstStab.verticies, cong_map, zero_odo}, width, height, false, debug);
    }

    firstOdo = firstStab.odometer;
    secondOdo = secondStab.odometer;

    return secondStab.verticies;
};

void gen(const int width, const int height, const int starting_width = 0, const int starting_height = 0, const bool debug = false) {
    int target_width, target_height;
    int SCALING = 2;
    std::vector<int> firstOdo;
    std::vector<int> secondOdo;
    std::vector<int> b;
    vector2D a;


    if (width == 1 || height == 1) {
        a = generateTorusCongruenceMap(width, height);
        b = generateNeutralElement(a, width, height, false, firstOdo, secondOdo, debug);
        createBMP(unflatten(b, width, height), nColourScheme(4), "../TorusOffset/"+std::to_string(width)+"x"+std::to_string(height)+".bmp");
        return;
    }

    if (width % 2 == 1) {
        target_width = width - 1;
    } else {
        target_width = width;
    }
    if (height % 2 == 1) {
        target_height = height - 1;
    } else {
        target_height = height;
    }

    int lowest_width = target_width;
    int t_width = lowest_width;
    int lowest_height = target_height;
    int t_height = lowest_height;

    while (t_width % SCALING == 0 && t_height % SCALING == 0) {
        t_width = t_width / SCALING;
        t_height = t_height / SCALING;
        if ( t_width >= starting_width && t_height >= starting_height) {
            lowest_width = t_width;
            lowest_height = t_height;
        }
    }

    if (starting_height) {
        firstOdo = readTXT("../TorusOffset/Odometers/"+std::to_string(starting_width)+"x"+std::to_string(starting_height)+"odo1.txt");
        secondOdo = readTXT("../TorusOffset/Odometers/"+std::to_string(starting_width)+"x"+std::to_string(starting_width)+"odo2.txt");
    } else {
        auto a = generateTorusCongruenceMap(lowest_width, lowest_height, 0);
        auto b = generateNeutralElement(a, lowest_width, lowest_height, false, firstOdo, secondOdo, debug);
    }

    for (int i = SCALING; i*lowest_width <= target_width; i = i*SCALING) {
        if (debug) { printf("\nGenerating %d x %d...\n", lowest_width * i, lowest_height * i);}
        firstOdo = scaleOdometer(firstOdo, lowest_width * (i/2), lowest_height * (i/2), SCALING);
        secondOdo = scaleOdometer(secondOdo, lowest_width * (i/2), lowest_height * (i/2), SCALING);
        a = generateTorusCongruenceMap(lowest_width * i, lowest_height * i, 0);
        b = generateNeutralElement(a, lowest_width * i, lowest_height * i, true, firstOdo, secondOdo, debug);
    }

    if (target_height != height || target_width != width){
        if (debug) { printf("\nGenerating %d x %d...\n", width, height);}
        a = generateTorusCongruenceMap(width, height, 0);
        firstOdo = increaseOdometer(firstOdo, target_width, target_height, 1);
        secondOdo = increaseOdometer(secondOdo, target_width, target_height, 1);
        b = generateNeutralElement( a, width, height, true, firstOdo, secondOdo, debug);
    }

    std::string filename1 = "../TorusOffset/Odometers/"+std::to_string(width)+"x"+std::to_string(height)+"odo1.txt";
    std::string filename2 = "../TorusOffset/Odometers/"+std::to_string(width)+"x"+std::to_string(height)+"odo2.txt";
                                
    createTXT(firstOdo, filename1);
    createTXT(secondOdo, filename2);
                            
    createBMP(unflatten(b, width, height), nColourScheme(4), "../TorusOffset/"+std::to_string(width)+"x"+std::to_string(height)+".bmp");

}

/* void plotOdo() {
    std::vector<int> firstOdo, secondOdo;
    int size = 400;


    auto b = doublePath(size, size, 2, true);

    firstOdo = firstOdo + secondOdo;

    // scale that shit!!!!

    float scaling = 255 / (float)(*std::max_element(firstOdo.begin(), firstOdo.end()) + 1);
    std::cout << scaling << "\n";

    for (int &i : firstOdo) {
        i = i * scaling;
    }
    
    createBMP(unflatten(firstOdo, size, size), 
                        nColourScheme(256), std::to_string(size)+"x"+std::to_string(size)+"odo.bmp");

}; */

int main(){
    srand(time(NULL));

    return 0;
};