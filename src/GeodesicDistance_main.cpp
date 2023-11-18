/*
  GeodesicDistance_main.cpp - Demonstrates the use of fp64Lib when calculating
  high precision floating point arithmetic with both Vincenty and Haversine
  geodesic distance formulas. Sample code provides timing to compare the
  precision and time cost of each function.

  @see https://fp64lib.org/ for documentation and details
  @see https://geodesyapps.ga.gov.au/vincenty-inverse for comparing calculated distances

*/

#include <Arduino.h>
#include <math.h>
#include <fp64lib.h>

// #############################################################################
// # GLOBAL CONSTANTS
// #############################################################################

/* (float64_t) Tampa, FL USA and Eiffel Tower, France */
const float64_t fp64_lat1 = fp64_sd(27.951692274858896);
const float64_t fp64_lon1 = fp64_sd(-82.45874393775438);
const float64_t fp64_lat2 = fp64_sd(48.85839511997923);
const float64_t fp64_lon2 = fp64_sd(2.294336581236128);

/* (double) Tampa, FL USA and Eiffel Tower, France */
const double lat1 = 27.951692274858896;
const double lon1 = -82.45874393775438;
const double lat2 = 48.85839511997923;
const double lon2 = 2.294336581236128;

// #############################################################################
// # FUNCTION DECLARATIONS
// #############################################################################

float64_t toRadians(float64_t deg);

double toRadians(double degrees);

float64_t getVincentyInverseDistance(float64_t lat1, float64_t lon1, float64_t lat2, float64_t lon2);

double getVincentyInverseDistance(double lat1, double lon1, double lat2, double lon2);

float64_t getHaversineDistance(float64_t lat1, float64_t lon1, float64_t lat2, float64_t lon2);

double getHaversineDistance(double lat1, double lon1, double lat2, double lon2);

// #############################################################################
// # ARDUINO SETUP
// #############################################################################

void setup()
{

  // Enable serial device
  Serial.begin(9600);
}

// #############################################################################
// # ARDUINO MAIN LOOP
// #############################################################################

void loop()
{

  long startTime = micros(), elapsedTime;

  float64_t vincentyDistanceFP64 = getVincentyInverseDistance(fp64_lat1, fp64_lon1, fp64_lat2, fp64_lon2);
  elapsedTime = micros() - startTime;
  Serial.println("Vincenty Distance (float64_t): " + String(fp64_to_string(fp64_div(vincentyDistanceFP64, fp64_sd(1000)), 18, 18)) + "(km)  Time to complete: " + elapsedTime + " us");
  Serial.println("Vincenty Distance (float64_t): " + String(fp64_to_string(vincentyDistanceFP64, 18, 18)) + "(m)  Time to complete: " + elapsedTime + " us\n");

  startTime = micros();
  double vincentyDistance = getVincentyInverseDistance(lat1, lon1, lat2, lon2);
  elapsedTime = micros() - startTime;
  Serial.println("Vincenty Distance    (double): " + String(vincentyDistance / 1000) + "(km)  Time to complete: " + elapsedTime + " us");
  Serial.println("Vincenty Distance    (double): " + String(vincentyDistance) + "(m)  Time to complete: " + elapsedTime + " us\n");

  startTime = micros();
  float64_t haversineDistanceFP64 = getHaversineDistance(fp64_lat1, fp64_lon1, fp64_lat2, fp64_lon2);
  elapsedTime = micros() - startTime;
  Serial.println("Haversine Distance (float64_t): " + String(fp64_to_string(fp64_div(haversineDistanceFP64, fp64_sd(1000)), 18, 18)) + "(km)  Time to complete: " + elapsedTime + " us");
  Serial.println("Haversine Distance (float64_t): " + String(fp64_to_string(haversineDistanceFP64, 18, 18)) + "(m)  Time to complete: " + elapsedTime + " us\n");

  startTime = micros();
  double haversineDistance = getHaversineDistance(lat1, lon1, lat2, lon2);
  elapsedTime = micros() - startTime;
  Serial.println("Haversine Distance   (double): " + String(haversineDistance / 1000) + "(km)  Time to complete: " + elapsedTime + " us");
  Serial.println("Haversine Distance   (double): " + String(haversineDistance) + "(m)  Time to complete: " + elapsedTime + " us\n\n\n");

  delay(60000);
}

// #############################################################################
// # FUNCTION IMPLEMENTATIONS
// #############################################################################

float64_t toRadians(float64_t deg)
{
  return fp64_mul(deg, fp64_div(float64_NUMBER_PI, fp64_sd(180.0)));
}

double toRadians(double degrees)
{
  return degrees * (PI / 180.0);
}

/*Calulcates the Vincenty inverse distance between two points using fp64lib */
float64_t getVincentyInverseDistance(float64_t lat1, float64_t lon1, float64_t lat2, float64_t lon2)
{

  lat1 = toRadians(lat1);
  lon1 = toRadians(lon1);
  lat2 = toRadians(lat2);
  lon2 = toRadians(lon2);

  float64_t a = fp64_sd(6378137);
  float64_t b = fp64_sd(6356752.314245);
  float64_t f = fp64_div(fp64_sd(1), fp64_sd(298.257223563));

  float64_t L = fp64_sub(lon2, lon1);

  float64_t U1 = fp64_atan(fp64_mul(fp64_sub(fp64_sd(1), f), fp64_tan(lat1)));
  float64_t U2 = fp64_atan(fp64_mul(fp64_sub(fp64_sd(1), f), fp64_tan(lat2)));
  float64_t sinU1 = fp64_sin(U1), cosU1 = fp64_cos(U1);
  float64_t sinU2 = fp64_sin(U2), cosU2 = fp64_cos(U2);

  float64_t lambda = L, lambdaP;
  float64_t sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, cos2SigmaM;

  int iterLimit = 100;

  do
  {
    float64_t sinLambda = fp64_sin(lambda), cosLambda = fp64_cos(lambda);
    sinSigma = fp64_sqrt(

        fp64_add(fp64_mul(fp64_mul(cosU2, sinLambda), fp64_mul(cosU2, sinLambda)),
                 fp64_mul((fp64_sub(fp64_mul(cosU1, sinU2), fp64_mul(fp64_mul(sinU1, cosU2), cosLambda))), (fp64_sub(fp64_mul(cosU1, sinU2), fp64_mul(fp64_mul(sinU1, cosU2), cosLambda))))));

    if (sinSigma == fp64_sd(0))
      return 0; // co-incident points

    cosSigma = fp64_add(fp64_mul(sinU1, sinU2), fp64_mul(fp64_mul(cosU1, cosU2), cosLambda));
    sigma = fp64_atan2(sinSigma, cosSigma);
    sinAlpha = fp64_div(fp64_mul(fp64_mul(cosU1, cosU2), sinLambda), sinSigma);
    cosSqAlpha = fp64_sub(fp64_sd(1), fp64_mul(sinAlpha, sinAlpha));
    cos2SigmaM = fp64_sub(cosSigma, fp64_div(fp64_mul(fp64_sd(2), fp64_mul(sinU1, sinU2)), cosSqAlpha));

    if (fp64_isnan(cos2SigmaM))
      cos2SigmaM = fp64_sd(0); // equatorial line: cosSqAlpha=0

    float64_t C = fp64_mul(fp64_mul(fp64_div(f, fp64_sd(16)), cosSqAlpha), (fp64_add(fp64_sd(4), fp64_mul(f, (fp64_sub(fp64_sd(4), fp64_mul(fp64_sd(3), cosSqAlpha)))))));
    lambdaP = lambda;
    lambda = fp64_add(L, fp64_mul((fp64_sub(fp64_sd(1), C)), fp64_mul(f, fp64_mul(sinAlpha,
                                                                                  (fp64_add(sigma, fp64_mul(C, fp64_mul(sinSigma, (fp64_add(cos2SigmaM, fp64_mul(C, fp64_mul(cosSigma,
                                                                                                                                                                             (fp64_add(fp64_sd(-1), fp64_mul(fp64_sd(2), fp64_mul(cos2SigmaM, cos2SigmaM))))))))))))))));
    int isGreater = fp64_compare(fp64_abs(fp64_sub(lambda, lambdaP)), fp64_sd(1e-12));

  } while (fp64_compare(fp64_abs(fp64_sub(lambda, lambdaP)), fp64_sd(1e-12)) > 0 && --iterLimit > 0);

  if (iterLimit == 0)
    return NAN; // formula failed to converge

  float64_t uSq = fp64_mul(cosSqAlpha, (fp64_div((fp64_sub(fp64_mul(a, a), fp64_mul(b, b))), (fp64_mul(b, b)))));
  float64_t A = fp64_add(fp64_sd(1), fp64_mul(fp64_div(uSq, fp64_sd(16384)), (fp64_add(fp64_sd(4096), fp64_mul(uSq, (fp64_add(fp64_sd(-768), fp64_mul(uSq, (fp64_sub(fp64_sd(320), fp64_mul(175, uSq)))))))))));
  float64_t B = fp64_mul(fp64_div(uSq, fp64_sd(1024)), (fp64_add(fp64_sd(256), fp64_mul(uSq, (fp64_add(fp64_sd(-128), fp64_mul(uSq, (fp64_sub(fp64_sd(74), fp64_mul(fp64_sd(47), uSq))))))))));
  float64_t deltaSigma = fp64_mul(B, fp64_mul(sinSigma, (fp64_add(cos2SigmaM, fp64_mul(fp64_div(B, fp64_sd(4)), (fp64_mul(cosSigma, (fp64_add(fp64_sd(-1), fp64_mul(fp64_sd(2), fp64_mul(cos2SigmaM, cos2SigmaM))) - fp64_mul(fp64_div(B, fp64_sd(6)), fp64_mul(cos2SigmaM, fp64_mul((fp64_add(fp64_sd(-3), fp64_mul(fp64_sd(4), fp64_mul(sinSigma, sinSigma)))), (fp64_add(fp64_sd(-3), fp64_mul(fp64_sd(4), fp64_mul(cos2SigmaM, cos2SigmaM)))))))))))))));

  float64_t s = fp64_mul(b, fp64_mul(A, (fp64_sub(sigma, deltaSigma))));

  return s;
}

/*Calulcates the Vincenty inverse distance between two points using double precision */
double getVincentyInverseDistance(double lat1, double lon1, double lat2, double lon2)
{

  lat1 = toRadians(lat1);
  lon1 = toRadians(lon1);
  lat2 = toRadians(lat2);
  lon2 = toRadians(lon2);

  double a = 6378137.0;
  double b = 6356752.314245;
  double f = 1 / 298.257223563;

  double L = lon2 - lon1;
  double U1 = atan((1 - f) * tan(lat1));
  double U2 = atan((1 - f) * tan(lat2));
  double sinU1 = sin(U1), cosU1 = cos(U1);
  double sinU2 = sin(U2), cosU2 = cos(U2);

  double lambda = L, lambdaP, iterLimit = 100.0;
  double cosSqAlpha, sinSigma, cos2SigmaM, cosSigma, sigma;
  double uSq, A, B, deltaSigma;
  do
  {
    double sinLambda = sin(lambda), cosLambda = cos(lambda);
    sinSigma = sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) +
                    (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
    if (sinSigma == 0)
      return 0; // co-incident points
    cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
    sigma = atan2(sinSigma, cosSigma);
    double sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
    cosSqAlpha = 1 - sinAlpha * sinAlpha;
    cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
    if (isnan(cos2SigmaM))
      cos2SigmaM = 0; // equatorial line: cosSqAlpha=0
    uSq = cosSqAlpha * (a * a - b * b) / (b * b);
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
    deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
    lambdaP = lambda;
    lambda = L + (1 - B) * f * sinAlpha * (sigma + B * sinSigma * (cos2SigmaM + B * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
  } while (abs(lambda - lambdaP) > 1e-12 && --iterLimit > 0);

  if (iterLimit == 0)
    return NAN; // formula failed to converge

  double s = b * A * (sigma - deltaSigma);

  return s;
}

/*Calulcates the Great Circle distance between two points using fp64lib */
float64_t getHaversineDistance(float64_t lat1, float64_t lon1, float64_t lat2, float64_t lon2)
{
  lat1 = toRadians(lat1);
  lon1 = toRadians(lon1);
  lat2 = toRadians(lat2);
  lon2 = toRadians(lon2);

  float64_t radiusOfEarth = fp64_sd(6378137.0); // Earth's radius in meters

  float64_t deltaLat = fp64_sub(lat2, lat1);
  float64_t deltaLon = fp64_sub(lon2, lon1);

  float64_t a = fp64_add(fp64_mul(fp64_sin(fp64_div(deltaLat, fp64_sd(2))), fp64_sin(fp64_div(deltaLat, fp64_sd(2)))), fp64_mul(fp64_cos(lat1), fp64_mul(fp64_cos(lat2), fp64_mul(fp64_sin(fp64_div(deltaLon, fp64_sd(2))), fp64_sin(fp64_div(deltaLon, fp64_sd(2)))))));

  float64_t c = fp64_mul(fp64_sd(2), fp64_atan2(fp64_sqrt(a), fp64_sqrt(fp64_sub(fp64_sd(1), a))));

  float64_t distance = fp64_mul(radiusOfEarth, c);

  return distance;
}

/*Calulcates the Great Circle distance between two points using double precision.*/
double getHaversineDistance(double lat1, double lon1, double lat2, double lon2)
{
  lat1 = toRadians(lat1);
  lon1 = toRadians(lon1);
  lat2 = toRadians(lat2);
  lon2 = toRadians(lon2);

  double radiusOfEarth = 6378137.0; // Earth's radius in meters

  double deltaLat = lat2 - lat1;
  double deltaLon = lon2 - lon1;

  double a = sin(deltaLat / 2) * sin(deltaLat / 2) +
             cos(lat1) * cos(lat2) *
                 sin(deltaLon / 2) * sin(deltaLon / 2);

  double c = 2 * atan2(sqrt(a), sqrt(1 - a));

  double distance = radiusOfEarth * c;

  return distance;
}
