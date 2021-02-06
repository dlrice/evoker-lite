from numpy import array, arctan2, pi, cos, sin, linspace, sqrt, dot
from numpy.linalg import eig
from scipy.stats import chi2
from struct import unpack
import logging


class SNPPosterior(object):
    """
    Based on tutorial at:
        http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
    """

    def __init__(self, file_path, n_batches):
        self.n_batches = n_batches
        self.open_file(file_path)
        self.BYTES_PER_BATCH = 4 * 33
        self.CALLS = ('11', '01', '00')  # = ('BB', 'AB', 'AA')

    def open_file(self, file_path):
        logging.debug(file_path)
        self.f = open(file_path, 'rb')

    def get_offset(self, variant_index, batch_index):
        """
        Affymetrix provide 33 parameters per SNP to describe the properties of the genotype calling
        ellipses for a batch (3,498 per SNP). Each parameter is represented by a 4-byte float. The 
        set of values for each SNP are ordered consecutively by batch (as given in the batch file).
        This is analagous to a matrix with rows=SNPs and columns=batch_parameters. The order of the 
        SNPs and batches are given by the BIM and Batch files.
        """
        offset = self.BYTES_PER_BATCH * \
            (variant_index * self.n_batches + batch_index)
        return offset

    def get_ellipses_parameters(self, variant_index, batch_index):
        """
        The 33 parameters are arranged in four groups. There are 7 parameters for each genotype cluster:
        BB, AB, AA and a fourth set of 12 parameters (CV) from Affymetrix's algorithm which can be ignored 
        for most analyses.

        The 7 parameters per cluster are five parameters (mu,sigma) to completely describe the ellipse and
        two counts of observations (NObsVar can be ignored for most analyses):

            1. mu_X
            2. sigma_00
            3. NObsMean
            4. NObsVar
            5. mu_Y
            6. sigma_11
            7. sigma_off-diagonal(01,10)

        Missing values are represented using the IEEE 754 NaN.
        """
        offset = self.get_offset(variant_index, batch_index)
        self.f.seek(offset)
        _bytes = self.f.read(self.BYTES_PER_BATCH)
        logging.debug('Reading: {}'.format(int(self.BYTES_PER_BATCH/4)))
        # 'f' is IEEE 754 binary32 or 4 bytes (=4*8=32)
        X = unpack('f'*int(self.BYTES_PER_BATCH/4), _bytes)
        X = array(X, dtype='float')
        # We only need 21 of the parameters as the spec states the other 12 can be ignored
        X = X[:21].reshape(3, -1)
        P = {}
        for parameters, call in zip(X, self.CALLS):
            mu_x, sigma_00, _, _, mu_y, sigma_11, sigma_01_10 = parameters
            P[call] = {
                'mean': array([mu_x, mu_y]),
                'covariance': array([[sigma_00, sigma_01_10], [sigma_01_10, sigma_11]]),
            }

        return P

    def get_batch_ellipse_points(self, variant_index, batch_index):
        parameters = self.get_ellipses_parameters(variant_index, batch_index)
        return {call: self.get_ellipse_points(**ps) for call, ps in parameters.items()}

    def get_ellipse_points(self, mean, covariance, confidence=0.85):
        eigenvalues, eigenvectors = eig(covariance)

        # Get the index of the largest eigenvector
        largest_eigenvector_index = eigenvalues.argmax()
        largest_eigenvalue = eigenvalues[largest_eigenvector_index]
        largest_eigenvector = eigenvectors[:, largest_eigenvector_index]

        # Get the smallest eigenvector and eigenvalue
        smallest_eigenvector_index = (largest_eigenvector_index + 1) % 2
        smallest_eigenvalue = eigenvalues[smallest_eigenvector_index]
        # smallest_eigenvector = eigenvectors[:, smallest_eigenvector_index]

        # Calculate the angle between the x-axis and the largest eigenvector
        angle = arctan2(largest_eigenvector[1], largest_eigenvector[0])

        # This angle is between -pi and pi.
        # Let's shift it such that the angle is between 0 and 2pi
        if angle < 0:
            angle += 2*pi

        # Get the confidence interval error ellipse
        chisquare_val = chi2.ppf(confidence, df=2)**(0.5)
        theta_grid = linspace(0, 2*pi)
        phi = angle
        X0 = mean[0]
        Y0 = mean[1]
        a = chisquare_val*sqrt(largest_eigenvalue)
        b = chisquare_val*sqrt(smallest_eigenvalue)

        # the ellipse in x and y coordinates
        ellipse_x_r = a*cos(theta_grid)
        ellipse_y_r = b*sin(theta_grid)

        # Define a rotation matrix
        R = array([[cos(phi), sin(phi)], [-sin(phi), cos(phi)]])

        # let's rotate the ellipse to some angle phi
        r_ellipse = dot(array([ellipse_x_r, ellipse_y_r]).T, R)
        return r_ellipse + array([X0, Y0])
