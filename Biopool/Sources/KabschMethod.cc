#include <KabschMethod.h>
//#include "Eigen/Geometry"

using namespace Victor;
using namespace Victor::Biopool;

// CONSTRUCTORS:

KabschMethod::KabschMethod() {
}

double KabschMethod::rotate(Spacer* set1, Spacer* set2) {
    Eigen::Affine3d output;
    //Matrice di rotazione
    output.linear() = Eigen::Matrix3d::Identity(3, 3);
    //Vettore di traslazione
    output.translation() = Eigen::Vector3d::Zero();
    int NumAmino1 = set1->sizeAmino();
    int NumAmino2 = set2->sizeAmino();
    Atom CAAtoms1[NumAmino1];
    Atom CAAtoms2[NumAmino2];
    Eigen::Matrix3Xd set1Matrix(3, NumAmino1);
    Eigen::Matrix3Xd set2Matrix(3, NumAmino2);

    for (int i = 0; i < NumAmino1; i++) {
        CAAtoms1[i] = set1->getAmino(i)[CA];
    }

    for (int i = 0; i < NumAmino2; i++) {
        CAAtoms2[i] = set2->getAmino(i)[CA];
    }

    for (int col = 0; col < NumAmino1; col++) {
        vgVector3<double> coords = CAAtoms1[col].getCoords();
        set1Matrix(0, col) = coords[0];
        set1Matrix(1, col) = coords[1];
        set1Matrix(2, col) = coords[2];
    }

    for (int col = 0; col < NumAmino2; col++) {
        vgVector3<double> coords = CAAtoms2[col].getCoords();
        set2Matrix(0, col) = coords[0];
        set2Matrix(1, col) = coords[1];
        set2Matrix(2, col) = coords[2];
    }


    //    //Questi sono i valori di ritorno, che includono la matrice di rotazione e il vettore di traslazione
    //    Eigen::Affine3d output;
    //    //Matrice di rotazione
    //    output.linear() = Eigen::Matrix3d::Identity(3, 3);
    //    //Vettore di traslazione
    //    output.translation() = Eigen::Vector3d::Zero();

    //Controlla che i due set abbiano uguale lunghezza
    if (set1Matrix.cols() != set2Matrix.cols())
        throw "I due set devono presentare uguale lunghezza";


    // Qui calcola le medie delle distanze e le sottrae poi per normalizzare i valori
    // INIZIO
    // First find the scale, by finding the ratio of sums of some distances,
    // then bring the datasets to the same scale.
    // Sta sottrendo le colonne a due a due
    double dist_in = 0, dist_out = 0;
    for (int col = 0; col < set1Matrix.cols() - 1; col++) {
        dist_in += (set1Matrix.col(col + 1) - set1Matrix.col(col)).norm();
        dist_out += (set2Matrix.col(col + 1) - set2Matrix.col(col)).norm();
    }

    //Restituisce l'output cosi' come' se non vi e' nessuna rotazione o traslazione da compiere
    if (dist_in <= 0 || dist_out <= 0)
        return 0; //MODIFICATO

    double scale = dist_out / dist_in;
    set2Matrix /= scale;

    // Find the centroids then shift to the origin
    Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
    Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
    for (int col = 0; col < set1Matrix.cols(); col++) {
        in_ctr += set1Matrix.col(col);
        out_ctr += set2Matrix.col(col);
    }
    in_ctr /= set1Matrix.cols();
    out_ctr /= set2Matrix.cols();
    for (int col = 0; col < set1Matrix.cols(); col++) {
        set1Matrix.col(col) -= in_ctr;
        set2Matrix.col(col) -= out_ctr;
    }
    //FINE

    // Calcolo della scomposizone svd
    Eigen::MatrixXd Cov = set1Matrix * set2Matrix.transpose();
    // Compie un tipo di svd computazionalmente piu' efficente di quella teorica
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Trova il verso in cui rotare per ottenere una rotazione destrorsa
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0)
        d = 1.0;
    else
        d = -1.0;
    //Genera la matrice I necessaria per trovare la rotazione
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    //Cambia, se necessario, l'ultimo valore per ottenere una rotazione verso destra
    I(2, 2) = d;
    Eigen::Matrix3d rotationMatrix = svd.matrixV() * I * svd.matrixU().transpose();

    //Ritorna la matrice di rotazione e il vettore di traslazione trovati
    output.linear() = scale * rotationMatrix;
    //Il vettore si calcola facendo la differenza tra il primo set rotato e il secondo
    output.translation() = scale * (out_ctr - rotationMatrix * in_ctr);
    return 0;
}

