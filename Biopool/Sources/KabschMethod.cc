#include <KabschMethod.h>


namespace Victor {
    namespace Biopool {

        // CONSTRUCTORS:

        KabschMethod::KabschMethod() {
        }

        void KabschMethod::rotate() {


//            //Questi sono i valori di ritorno, che includono la matrice di rotazione e il vettore di traslazione
//            Eigen::Affine3d output;
//            //Matrice di rotazione
//            output.linear() = Eigen::Matrix3d::Identity(3, 3);
//            //Vettore di traslazione
//            output.translation() = Eigen::Vector3d::Zero();
//
//            //Controlla che i due set abbiano uguale lunghezza
//            if (set1.cols() != set2.cols())
//                throw "I due set devono presentare uguale lunghezza";
//
//
//            // Qui calcola le medie delle distanze e le sottrae poi per normalizzare i valori
//            // INIZIO
//            // First find the scale, by finding the ratio of sums of some distances,
//            // then bring the datasets to the same scale.
//            // Sta sottrendo le colonne a due a due
//            double dist_in = 0, dist_out = 0;
//            for (int col = 0; col < set1.cols() - 1; col++) {
//                dist_in += (set1.col(col + 1) - set1.col(col)).norm();
//                dist_out += (set2.col(col + 1) - set2.col(col)).norm();
//            }
//
//            //Restituisce l'output cosi' come' se non vi e' nessuna rotazione o traslazione da compiere
//            if (dist_in <= 0 || dist_out <= 0)
//                return output;
//
//            double scale = dist_out / dist_in;
//            set2 /= scale;
//
//            // Find the centroids then shift to the origin
//            Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
//            Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
//            for (int col = 0; col < set1.cols(); col++) {
//                in_ctr += set1.col(col);
//                out_ctr += set2.col(col);
//            }
//            in_ctr /= set1.cols();
//            out_ctr /= set2.cols();
//            for (int col = 0; col < set1.cols(); col++) {
//                set1.col(col) -= in_ctr;
//                set2.col(col) -= out_ctr;
//            }
//            //FINE
//
//            // Calcolo della scomposizone svd
//            Eigen::MatrixXd Cov = set1 * set2.transpose();
//            // Compie un tipo di svd computazionalmente piu' efficente di quella teorica
//            Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
//
//            // Trova il verso in cui rotare per ottenere una rotazione destrorsa
//            double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
//            if (d > 0)
//                d = 1.0;
//            else
//                d = -1.0;
//            //Genera la matrice I necessaria per trovare la rotazione
//            Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
//            //Cambia, se necessario, l'ultimo valore per ottenere una rotazione verso destra
//            I(2, 2) = d;
//            Eigen::Matrix3d rotationMatrix = svd.matrixV() * I * svd.matrixU().transpose();
//
//            //Ritorna la matrice di rotazione e il vettore di traslazione trovati
//            output.linear() = scale * rotationMatrix;
//            //Il vettore si calcola facendo la differenza tra il primo set rotato e il secondo
//            output.translation() = scale * (out_ctr - rotationMatrix * in_ctr);
        }

    }
} // namespace
