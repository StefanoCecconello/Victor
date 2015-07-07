#include <KabschMethod.h>
#include "Eigen/Geometry"

using namespace Victor;
using namespace Victor::Biopool;

// CONSTRUCTORS:

KabschMethod::KabschMethod() {
}
/**
 * This method implement the KabschMethod method for obtein the rototraslation 
 * that minimize the rmsd value of the rotated protein and the second protein.
 * 
 * @param set1Matrix (Eigen::Matrix3Xd), the 3*N matrix with the coordinate
 * of the atom of the first protein;
 * @param set2Matrix (Eigen::Matrix3Xd), the 3*N matrix with the coordinate
 * of the atom of the second protein;
 * @return Eigen::Affine3d*, the rototraslation returned by the method.
 */
Eigen::Affine3d* KabschMethod::rotate(Eigen::Matrix3Xd set1Matrix, Eigen::Matrix3Xd set2Matrix) const{
    Eigen::Affine3d* output = new Eigen::Affine3d();
    //Rotation Matrix
    output->linear() = Eigen::Matrix3d::Identity(3, 3);
    //Translation Vector
    output->translation() = Eigen::Vector3d::Zero();



    //Questi sono i valori di ritorno, che includono la matrice di rotazione e il vettore di traslazione


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
        return output;

    //double scale = dist_out / dist_in;
    //set2Matrix /= scale;

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
    //output->linear() = scale * rotationMatrix;
    output->linear() = rotationMatrix;
    //Il vettore si calcola facendo la differenza tra il primo set rotato e il secondo

    //output->translation() = scale * (out_ctr - rotationMatrix * in_ctr);
    output->translation() = out_ctr - rotationMatrix * in_ctr;
    return output;
}

