#include "sparse_matrix.h"
#include "metrics.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include "debug.h"

using namespace std;

typedef vector<vector<unsigned char> > Image;

static const double epsilon = numeric_limits<double>::epsilon();

struct SimulationData
{
    Image image;
    unsigned cell_size;
    unsigned discr_size;
    unsigned method;
    vector<double> noise_levels;
};

bool load_csv_image(const string& filename, Image& image)
{
    ifstream ifile(filename);
    if (ifile.fail()) {
        return false;
    }

    string line;
    while (getline(ifile, line)) {
        image.push_back(vector<unsigned char>());
        stringstream ss(line);
        string number;
        while (getline(ss, number, ',')) {
            image.back().push_back(stoi(number));
        }
    }

    return true;
}

unsigned char convert_to_pixel(double x)
{
    if (x > 255.0) {
        x = 255.0;
    }
    else if (x < 0.0) {
        x = 0.0;
    }
    
    return (unsigned char)x;
}

vector<Image> convert_to_images(const vector<Vector>& s, unsigned rowsize)
{
    vector<Image> res(s.size(), Image(rowsize, vector<unsigned char>(rowsize)));
    for (unsigned k = 0; k < s.size(); k++) {
        for (unsigned i = 0; i < rowsize; i++) {
            for (unsigned j = 0; j < rowsize; j++) {
                res[k][i][j] = convert_to_pixel(s[k][i*rowsize + j]);
            }
        }
    }
    return res;
}

Image scale(const Image& img, unsigned size, unsigned cell_size)
{
    Image res(size, vector<unsigned char>(size));
    for (unsigned i = 0; i < size; i++) {
        for (unsigned j = 0; j < size; j++) {
            res[i][j] = img[i/cell_size][j/cell_size];
        }
    }
    return res;
}

double get_psnr(const Image& img1, const Image& img2)
{
    unsigned n = img1.size();
    double ecm = 0.0;
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++) {
            ecm += (img1[i][j] - img2[i][j])*(img1[i][j] - img2[i][j]);
        }
    }
    ecm /= n*n;
    
    return 10.0 * log10((255.0*255.0) / ecm);
}

void save_as_csv_image(const string& filename, const Image& img)
{
    ofstream ofile(filename);

    for (unsigned i = 0; i < img.size(); i++) {
        ofile << (unsigned)(img[i][0]);
        for (unsigned j = 1; j < img[i].size(); j++) {
            ofile << ", " << (unsigned)(img[i][j]);
        }
        ofile << endl;
    }
}

unsigned celda(unsigned x, unsigned y, const SimulationData& sd)
{
    return ((y/sd.cell_size) * sd.discr_size + x / sd.cell_size);
}

void simulate_ray
(
    const SimulationData& sd,
    unsigned ray_index,
    unsigned x0,
    unsigned y0,
    unsigned x1,
    unsigned y1,
    SparseMatrix& D,
    vector<Vector>& ts
)
{
    double time = 0.0;
    vector<double> distances(sd.discr_size*sd.discr_size);
    
    // Un rayo se representa con una funcion lineal y=ax+b, donde el eje y va de arriba hacia abajo
    // La funcion lineal debe pasar por los centros de los pixeles (x0,y0) y (x1,y1)
    
    // Coordenadas exactas de los puntos por donde pasa la funcion:
    double x0d = x0 + 0.5;
    double y0d = y0 + 0.5;
    double x1d = x1 + 0.5;
    double y1d = y1 + 0.5;
    
    unsigned posx;
    unsigned posy;
    if (x0 < x1) {
        posx = x0;
        posy = y0;
    }
    else if (x1 < x0) {
        posx = x1;
        posy = y1;
    }
    else {
        // Si el rayo es vertical, lo "torcemos" un poco
        // Por convencion empezamos en el pixel (x0,y0)
        posx = x0;
        posy = y0;
        x0d -= 0.25;
        x1d += 0.25;
    }
    
    // Determinamos los coeficientes de la funcion
    double a = (y1d - y0d)/(x1d - x0d);
    double b = y0d - a*x0d;
    
    
    while (posx < sd.image.size() && posy < sd.image.size()) { //mientras no me salga de la imagen

        double actual = (double)(sd.image[posy][posx]);

        //sumo actual + 1 al tiempo total de este rayo k
        time += actual;

        //ahora sumo uno a la distancia en esta celda
        distances[celda(posx, posy, sd)] += 1.0;

        double derecha = a*((double)(posx + 1)) + b; //evaluo en el borde derecho del pixel actual

        if (fabs(derecha - (double)posy) < epsilon) {
            posy--;
            posx++;
        }
        else if (fabs(derecha - (double)(posy+1)) < epsilon) {
            posx++;
            posy++;
        }
        else if (derecha < (double)posy) {
            posy--;
        }
        else if ((double)posy < derecha && derecha < (double)(posy + 1)) {
            posx++;
        }
        else { // if ((double)(posy + 1) < derecha)
            posy++;
        }
    }
    
    for (unsigned j = 0; j < D.num_columns(); j++) {
        if (distances[j] != 0.0) {
            D.get_column(j).push_back(make_pair(ray_index, distances[j]));
        }
    }

    // Guardamos los tiempos tardados (agregando ruido)
    for (unsigned i = 0; i < ts.size(); i++) {
        double noise = ((2.0*sd.noise_levels[i])/RAND_MAX)*rand() - sd.noise_levels[i];
        ts[i][ray_index] = (time + noise) >= 0.0 ? time + noise : 0.0;
    }
}

enum Direction {
    UP_LEFT,
    UP_RIGHT,
    DOWN_LEFT,
    DOWN_RIGHT
};

/* Esto es para generar un rayo diagonal de pendiente 1 sin 
 * complicarse la vida.
 * Hace "trampa" porque llama a simulate_ray con un pixel que 
 * no va a estar en el borde de la imagen, pero funciona. */
void simulate_ray
(
    const SimulationData& sd,
    unsigned ray_index,
    unsigned x0,
    unsigned y0,
    Direction dir,
    SparseMatrix& D,
    vector<Vector>& ts
)
{
    switch (dir) {
    case UP_LEFT:
        simulate_ray(sd, ray_index, x0, y0, x0 - 1, y0 - 1, D, ts);
        break;
    case UP_RIGHT:
        simulate_ray(sd, ray_index, x0, y0, x0 + 1, y0 - 1, D, ts);
        break;
    case DOWN_LEFT:
        simulate_ray(sd, ray_index, x0, y0, x0 - 1, y0 + 1, D, ts);
        break;
    case DOWN_RIGHT:
        simulate_ray(sd, ray_index, x0, y0, x0 + 1, y0 + 1, D, ts);
        break;
    default:
        break;
    }
}

void simulate(const SimulationData& sd, SparseMatrix& D, vector<Vector>& ts)
{
    cout << "Simulando tomografia..." << endl;
    
    unsigned num_cells = sd.discr_size * sd.discr_size;
    unsigned imgsize = sd.image.size();
    
    // Generar todos los posibles rayos que partan de un lado y lleguen al lado opuesto
    if (sd.method == 0) {
        unsigned num_rays = 2 * imgsize * imgsize;
        D = SparseMatrix(num_rays, num_cells);
        
        for (unsigned i = 0; i < ts.size(); i++) {
            ts[i].resize(num_rays);
        }
        
        unsigned ray_index = 0;
        
        for (unsigned y0 = 0; y0 < imgsize; y0++) {
            for (unsigned y1 = 0; y1 < imgsize; y1++) {
                simulate_ray(sd, ray_index, 0, y0, imgsize - 1, y1, D, ts);
                ray_index++;
            }
        }
        
        for (unsigned x0 = 0; x0 < imgsize; x0++) {
            for (unsigned x1 = 0; x1 < imgsize; x1++) {
                simulate_ray(sd, ray_index, x0, 0, x1, imgsize - 1, D, ts);
                ray_index++;
            }
        }
    }
    
    // Rayos verticales, horizontales y diagonales
    else if (sd.method == 1) {
        unsigned num_rays = 6 * imgsize - 6;
        D = SparseMatrix(num_rays, num_cells);
        
        for (unsigned i = 0; i < ts.size(); i++) {
            ts[i].resize(num_rays);
        }
        
        unsigned ray_index = 0;
        
        for (unsigned y = 0; y < imgsize; y++) {
            simulate_ray(sd, ray_index, 0, y, imgsize - 1, y, D, ts);
            ray_index++;
        }
        for (unsigned x = 0; x < imgsize; x++) {
            simulate_ray(sd, ray_index, x, 0, x, imgsize - 1, D, ts);
            ray_index++;
        }
        
        for (unsigned y = 0; y < imgsize - 1; y++) {
            simulate_ray(sd, ray_index, 0, y, DOWN_RIGHT, D, ts);
            ray_index++;
        }
        
        for (unsigned x = 1; x < imgsize - 1; x++) {
            simulate_ray(sd, ray_index, x, 0, DOWN_RIGHT, D, ts);
            ray_index++;
        }
        
        for (unsigned y = 1; y < imgsize; y++) {
            simulate_ray(sd, ray_index, 0, y, UP_RIGHT, D, ts);
            ray_index++;
        }
        
        for (unsigned x = 1; x < imgsize - 1; x++) {
            simulate_ray(sd, ray_index, x, imgsize - 1, UP_RIGHT, D, ts);
            ray_index++;
        }
    }
    
    // Desde cada una de las cuatro esquinas barrer toda la imagen con rayos
    else if (sd.method == 2) {
        unsigned num_rays = 4*(2 * imgsize - 1);
        D = SparseMatrix(num_rays, num_cells);
        
        for (unsigned i = 0; i < ts.size(); i++) {
            ts[i].resize(num_rays);
        }
        
        unsigned ray_index = 0;
        
        for (unsigned y = 0; y < imgsize; y++) {
            simulate_ray(sd, ray_index, 0, 0, imgsize - 1, y, D, ts);
            ray_index++;
        }
        for (unsigned x = 0; x < imgsize - 1; x++) {
            simulate_ray(sd, ray_index, 0, 0, x, imgsize - 1, D, ts);
            ray_index++;
        }
        
        for (unsigned y = 0; y < imgsize; y++) {
            simulate_ray(sd, ray_index, imgsize - 1, 0, 0, y, D, ts);
            ray_index++;
        }
        for (unsigned x = 1; x < imgsize; x++) {
            simulate_ray(sd, ray_index, imgsize - 1, 0, x, imgsize - 1, D, ts);
            ray_index++;
        }
        
        for (unsigned y = 0; y < imgsize; y++) {
            simulate_ray(sd, ray_index, 0, imgsize - 1, imgsize - 1, y, D, ts);
            ray_index++;
        }
        for (unsigned x = 0; x < imgsize - 1; x++) {
            simulate_ray(sd, ray_index, 0, imgsize - 1, x, 0, D, ts);
            ray_index++;
        }
        
        for (unsigned y = 0; y < imgsize; y++) {
            simulate_ray(sd, ray_index, imgsize - 1, imgsize - 1, 0, y, D, ts);
            ray_index++;
        }
        for (unsigned x = 1; x < imgsize; x++) {
            simulate_ray(sd, ray_index, imgsize - 1, imgsize - 1, x, 0, D, ts);
            ray_index++;
        }
    }
    
    // Rayos aleatorios, sd.method indica la cantidad de rayos a generar
    else {
      
        unsigned num_rays = sd.method;
        D = SparseMatrix(num_rays, num_cells);
        
        for (unsigned i = 0; i < ts.size(); i++) {
            ts[i].resize(num_rays);
        }
        
        for (unsigned i = 0; i < num_rays; i++) {
            unsigned r0 = rand() % 6;
            unsigned r1 = rand() % imgsize;
            unsigned r2 = rand() % imgsize;
            if (r0 == 0) {
                simulate_ray(sd, i, 0, r1, imgsize - 1, r2, D, ts);
            }
            else if (r0 == 1) {
                simulate_ray(sd, i, r1, 0, r2, imgsize - 1, D, ts);
            }
            else if (r0 == 2) {
                simulate_ray(sd, i, 0, r1, r2, 0, D, ts);
            }
            else if (r0 == 3) {
                simulate_ray(sd, i, r1, 0, imgsize - 1, r2, D, ts);
            }
            else if (r0 == 4) {
                simulate_ray(sd, i, imgsize - 1, r1, r2, imgsize - 1, D, ts);
            }
            else {
                simulate_ray(sd, i, r1, imgsize - 1, 0, r2, D, ts);
            }
        }
    }
}

void output_results(const SimulationData& sd, const Metrics& metrics)
{
    ofstream ofile("results.txt", std::ios::app);
    ofile << sd.image.size() << " " << sd.cell_size << " " << sd.discr_size << " " << sd.method << endl;
    ofile << metrics.cond_number << endl;
    ofile << metrics.reconstruction_time << endl;
    ofile << metrics.num_eigen_found << endl;
    ofile << sd.noise_levels.size() << endl;
    for (unsigned i = 0; i < sd.noise_levels.size(); i++) {
        ofile << sd.noise_levels[i] << " " << metrics.psnr[i] << endl;
    }
    
    ofile.close();
}

int main(int argc, char* argv[])
{
    if (argc < 6) {
        cout << "Error: parametros invalidos." << endl;
        return 1;
    }
    
    srand(1000);
    
    SimulationData sd;
    Metrics metrics;
    
    // Leemos parametros
    string img_name_in = argv[1];
    string img_name_out = argv[2];
    sd.cell_size = stoi(argv[3]);
    sd.method = stoi(argv[4]);

    // Leemos los niveles de ruido y armamos los nombres de salida
    vector<string> out_names;
    int i = 5;
    while (i < argc) {
        sd.noise_levels.push_back(atof(argv[i]));
        out_names.push_back(img_name_out);
        size_t pos = out_names.back().rfind('.');
        if (pos != string::npos) {
            out_names.back().insert(pos, string("_") + argv[i]);
        }
        else {
            out_names.back().append(string("_") + argv[i]);
        }
        i++;
    }
    
    metrics.psnr.resize(sd.noise_levels.size());

    // Cargamos la imagen de entrada
    if (!load_csv_image(img_name_in, sd.image)) {
        cout << "Error: no se pudo abrir " << img_name_in << "." << endl;
        return 1;
    }
    
    // Obtenemos tamaño de la discretizacion
    sd.discr_size = sd.image.size() / sd.cell_size;
    if (sd.image.size() % sd.cell_size != 0) {
        sd.discr_size++;
    }

    // Simulamos la tomografia, obteniendo la matriz D y los vectores t (y el tiempo de ejecucion)
    SparseMatrix D;
    vector<Vector> ts(sd.noise_levels.size());
    simulate(sd, D, ts);
    
    // En este punto se puede imprimir la matriz D en un archivo, con 
    // la funcion print de debug.h, para luego ver los autovalores de DtD 
    // con el script print_eigenvalues.py
    
    // Reconstruimos la imagen
    cout << "Reconstruyendo imagen..." << endl;
    vector<Vector> s = least_squares(D, ts, metrics);
    vector<Image> results = convert_to_images(s, sd.discr_size);
    for (unsigned i = 0; i < results.size(); i++) {
        save_as_csv_image(out_names[i], results[i]);
        metrics.psnr[i] = get_psnr(scale(results[i], sd.image.size(), sd.cell_size), sd.image);
    }
    
    cout << "Tiempo de reconstruccion: " << metrics.reconstruction_time << " segundos." << endl;
    cout << "Numero de condicion de la matriz DtD: " << metrics.cond_number << endl;
    
    for (unsigned i = 0; i < metrics.psnr.size(); i++) {
        cout << "PSNR correspondiente a nivel de ruido " << sd.noise_levels[i] << ": " << metrics.psnr[i] << endl;
    }
    
    // output_results(sd, metrics); // Esto agrega los resultados obtenidos a un archivo de texto

    return 0;
}
