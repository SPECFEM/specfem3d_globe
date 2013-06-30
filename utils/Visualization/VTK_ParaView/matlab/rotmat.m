function  rotation_matrix = rotmat(clon,clat,grot)
% this function calculate the 3x3 rotation matrix from the AB chunk
% frame to the actual frame defined by (clon,clat,grot)

% compute colatitude and longitude and convert to radians
  DEGREES_TO_RADIANS = pi / 180.0;
  alpha = clon * DEGREES_TO_RADIANS;
  beta = (90.00 - clat) * DEGREES_TO_RADIANS;
  gamma = grot * DEGREES_TO_RADIANS;

  sina = sin(alpha);
  cosa = cos(alpha);
  sinb = sin(beta);
  cosb = cos(beta);
  sing = sin(gamma);
  cosg = cos(gamma);

% define rotation matrix
  rotation_matrix(1,1) = cosg*cosb*cosa-sing*sina;
  rotation_matrix(1,2) = -sing*cosb*cosa-cosg*sina;
  rotation_matrix(1,3) = sinb*cosa;
  rotation_matrix(2,1) = cosg*cosb*sina+sing*cosa;
  rotation_matrix(2,2) = -sing*cosb*sina+cosg*cosa;
  rotation_matrix(2,3) = sinb*sina;
  rotation_matrix(3,1) = -cosg*sinb;
  rotation_matrix(3,2) = sing*sinb;
  rotation_matrix(3,3) = cosb;

