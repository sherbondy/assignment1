#ifdef __APPLE__
#include <OpenGL/gl.h>
/* Just in case we need these later
// References:
// http://alumni.cs.ucsb.edu/~wombatty/tutorials/opengl_mac_osx.html
// # include <OpenGL/gl.h>
// # include <OpenGL/glu.h>
*/
#else
#include <GL/gl.h>
#endif

#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
using namespace std;

namespace
{
    // Approximately equal to.  We don't want to use == because of
    // precision issues with floating point.
    inline bool approx( const Vector3f& lhs, const Vector3f& rhs )
    {
        const float eps = 1e-8f;
        return ( lhs - rhs ).absSquared() < eps;
    }

    // The Bernstein matrix and its cuddly derivatives wrt time.
    
    const Matrix4f BEZ(1, -3,  3, -1,
                       0,  3, -6,  3,
                       0,  0,  3, -3,
                       0,  0,  0,  1);
    
    const Matrix4f DBEZ(-3,   6, -3, 0,
                         3, -12,  9, 0,
                         0,   6, -9, 0,
                         0,   0,  3, 0);
    
    // not actually used in this assignment
    const Matrix4f D2_BEZ(  6,   -6, 0, 0,
                          -12,   18, 0, 0,
                            6,  -18, 0, 0,
                            0,    6, 0, 0);
    
    const Matrix4f BEZ_INV(1,     1,     1, 1,
                           0, 1.f/3, 2.f/3, 1,
                           0,     0, 1.f/3, 1,
                           0,     0,     0, 1);
    
    const Matrix4f BSP(1.f/6, -1.f/2, 1.f/2, -1.f/6,
                       2.f/3,      0,    -1,  1.f/2,
                       1.f/6,  1.f/2, 1.f/2, -1.f/2,
                       0,          0,     0,  1.f/6);
    
    Matrix4f GforPAtIndex(const vector < Vector3f >& P, unsigned i){
        Vector4f P1(P[i],   0);
        Vector4f P2(P[i+1], 0);
        Vector4f P3(P[i+2], 0);
        Vector4f P4(P[i+3], 0);
        
        Matrix4f G = Matrix4f(P1, P2, P3, P4, true);
        return G;
    }
}

Curve evalBezier( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 || P.size() % 3 != 1 )
    {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit( 0 );
    }

    // TODO:
    // You should implement this function so that it returns a Curve
    // (e.g., a vector< CurvePoint >).  The variable "steps" tells you
    // the number of points to generate on each piece of the spline.
    // At least, that's how the sample solution is implemented and how
    // the SWP files are written.  But you are free to interpret this
    // variable however you want, so long as you can control the
    // "resolution" of the discretized spline curve with it.

    // Make sure that this function computes all the appropriate
    // Vector3fs for each CurvePoint: V,T,N,B.
    // [NBT] should be unit and orthogonal.

    // Also note that you may assume that all Bezier curves that you
    // receive have G1 continuity.  Otherwise, the TNB will not be
    // be defined at points where this does not hold.
    
    // Premise: create a bezier curve for each chunk of 4 control points
    // Then take steps samples from the curve by going at t increments of 1/steps
    // make sure to normalize tangents/normals/binormals
    
    Curve R;

    cerr << "\t>>> evalBezier has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        P[i].print();
    }
    
    for (unsigned i = 0; i < (P.size() - 1); i += 3) {
        // make a 4x4 matrix of 4 control points, plus a row of 0s.
        // when creating the vec3fs for the points, ignore the 4th coord.
        
        // setting the columns here
        // appending 0s to the point ends to let us make a 4x4 matrix
        Matrix4f G = GforPAtIndex(P, i);
        
        Matrix4f GB  = G * BEZ;  // to compute V
        Matrix4f GdB = G * DBEZ; // to compute T
                
        for (unsigned step = 0; step <= steps; step += 1) {
            float t = (1.0/steps) * float(step);
            Vector4f times(1, t, pow(t, 2), pow(t, 3));
            
            // We have to call .xyz() because we used a 4x4 matrix for
            // G instead of a 3x4
            Vector3f V = (GB  * times).xyz();
            Vector3f T = (GdB * times).xyz().normalized();
            
            Vector3f N, B;
            
            if (step == 0) {
                // arbitrarily set B to the forward vector
                B = Vector3f::FORWARD;
                if (Vector3f::cross(T, B) == Vector3f::ZERO) {
                    // change B if first selection was parallel to T.
                    B = Vector3f::UP;
                }
                
                N = Vector3f::cross(B, T).normalized();
            } else {
                CurvePoint previous = R.back();
                
                N = Vector3f::cross(previous.B, T).normalized();
                B = Vector3f::cross(T, N).normalized();
            }
            
            CurvePoint newPoint = { V, T, N, B};
            R.push_back(newPoint);
        }
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;

    return R;
}

Curve evalBspline( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 )
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit( 0 );
    }

    // TODO:
    // It is suggested that you implement this function by changing
    // basis from B-spline to Bezier.  That way, you can just call
    // your evalBezier function.
    
    // Q: how do I chunk the control points to make suitable beziers?
    // A: take 4, shift 1, take 4, shift 1, until your reach the end.
    
    // multiply groups of 4 points by BSP*BEZ_INV to yield equivalent
    // control points in bezier space.

    cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): " << P.size() << endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        P[i].print();
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;
    
    vector< Vector3f > transformedPoints;
    for (unsigned i = 0; i < (P.size() - 1); i += 3) {
        Matrix4f G = GforPAtIndex(P, i);
        
        Matrix4f GBezBspInv = G * BSP * BEZ_INV;
        for (unsigned col = 0; col < 4; ++col) {
            Vector3f transPoint = GBezBspInv.getCol(col).xyz();
            transformedPoints.push_back(transPoint);
        }
    }
    
    Curve bspR = evalBezier(transformedPoints, steps);

    // Return an empty curve right now.
    return bspR;
}

Curve evalCircle( float radius, unsigned steps )
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).
    
    // Preallocate a curve with steps+1 CurvePoints
    Curve R( steps+1 );

    // Fill it in counterclockwise
    for( unsigned i = 0; i <= steps; ++i )
    {
        // step from 0 to 2pi
        float t = 2.0f * M_PI * (float( i ) / steps);

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f( cos(t), sin(t), 0 );
        
        // Tangent vector is first derivative
        R[i].T = Vector3f( -sin(t), cos(t), 0 );
        
        // Normal vector is second derivative
        R[i].N = Vector3f( -cos(t), -sin(t), 0 );

        // Finally, binormal is facing up.
        R[i].B = Vector3f( 0, 0, 1 );
    }

    return R;
}

void drawCurve( const Curve& curve, float framesize )
{
    // Save current state of OpenGL
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    // Setup for line drawing
    glDisable( GL_LIGHTING ); 
    glColor4f( 1, 1, 1, 1 );
    glLineWidth( 1 );
    
    // Draw curve
    glBegin( GL_LINE_STRIP );
    for( unsigned i = 0; i < curve.size(); ++i )
    {
        glVertex( curve[ i ].V );
    }
    glEnd();

    glLineWidth( 1 );

    // Draw coordinate frames if framesize nonzero
    if( framesize != 0.0f )
    {
        Matrix4f M;

        for( unsigned i = 0; i < curve.size(); ++i )
        {
            M.setCol( 0, Vector4f( curve[i].N, 0 ) );
            M.setCol( 1, Vector4f( curve[i].B, 0 ) );
            M.setCol( 2, Vector4f( curve[i].T, 0 ) );
            M.setCol( 3, Vector4f( curve[i].V, 1 ) );

            glPushMatrix();
            glMultMatrixf( M );
            glScaled( framesize, framesize, framesize );
            glBegin( GL_LINES );
            glColor3f( 1, 0, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 1, 0, 0 );
            glColor3f( 0, 1, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 1, 0 );
            glColor3f( 0, 0, 1 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 0, 1 );
            glEnd();
            glPopMatrix();
        }
    }
    
    // Pop state
    glPopAttrib();
}

