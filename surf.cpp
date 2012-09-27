#include "surf.h"
#include "extra.h"
using namespace std;

namespace
{
    
    // We're only implenting swept surfaces where the profile curve is
    // flat on the xy-plane.  This is a check function.
    static bool checkFlat(const Curve &profile)
    {
        for (unsigned i=0; i<profile.size(); i++) {
            if (profile[i].V[2] != 0.0 ||
                profile[i].T[2] != 0.0 ||
                profile[i].N[2] != 0.0) {
                return false;
            }
        }
    
        return true;
    }
}

void addFacesToSurface(Surface& surface, unsigned i, unsigned step,
                       unsigned steps, unsigned profileSize)
{
    // do not create additional triangles once we're at the bottom of the curve
    if ( i + 1 != profileSize) {
        unsigned vIndex     = step*profileSize + i;
        // make sure next loops back to zero if you're on the last segment
        unsigned nextVIndex = (vIndex + profileSize) % (steps * profileSize);
        
        Tup3u face1(vIndex, vIndex + 1, nextVIndex);
        Tup3u face2(vIndex + 1, nextVIndex + 1, nextVIndex);
        surface.VF.push_back(face1);
        surface.VF.push_back(face2);
    }
}

void sweepProfile(const Curve &profile, unsigned profileSize, unsigned step,
                  unsigned steps, Surface &surface, Matrix4f transform)
{
    
    Matrix3f transformInvT = transform.getSubmatrix3x3(0, 0).inverse().transposed();
    
    for (unsigned i = 0; i < profileSize; ++i) {
        CurvePoint profilePoint = profile[i];
        Vector3f profileB = profilePoint.B;
        
        Vector3f vertexT = (transform * Vector4f(profilePoint.V, 1)).xyz();
        // now the inverse transpose bit actually matters for drawing normals
        Vector3f normalT = (transformInvT * profilePoint.N).xyz().normalized();
        
        surface.VV.push_back(vertexT);
        surface.VN.push_back(normalT);
        
        addFacesToSurface(surface, i, step, steps, profileSize);
    }
}

Surface makeSurfRev(const Curve &profile, unsigned steps)
{
    Surface surface;
    
    if (!checkFlat(profile))
    {
        cerr << "surfRev profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // idea: copy profile, transformed over and over by increment of rotation matrix.
    // faces: curve point 0, same curve point 1, adjacent curve point 0.
    //   and: curve point 1, adjacent curve point 1, adjacent curve point 0.
    // If I try making my own surfaces, always put profile on left of y-axis (make all x-coords -)
    
    // RotateY(radians) * CurvePoint.V, do for each point in profile, radians = 2pi/steps
    
    // first iteration: establish the vectors and normals
    unsigned profileSize = profile.size();
    
    for (unsigned step = 0; step < steps; step++) {
        float rotation = 2 * M_PI * step / steps;
        Matrix4f MRotY = Matrix4f::rotateY(rotation);
        
        sweepProfile(profile, profileSize, step, steps, surface, MRotY);
    }
     
    return surface;
}

Surface makeGenCyl(const Curve &profile, const Curve &sweep )
{
    Surface surface;

    if (!checkFlat(profile))
    {
        cerr << "genCyl profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    unsigned sweepSize   = sweep.size();
    unsigned profileSize = profile.size();

    for (unsigned step = 0; step < sweepSize; ++step) {
        CurvePoint point = sweep[step];
        // column-wise specification of the transformation
        Matrix4f transform(Vector4f(point.N, 0),
                           Vector4f(point.B, 0),
                           Vector4f(point.T, 0),
                           Vector4f(point.V, 1),
                           true);
                
        sweepProfile(profile, profileSize, step, sweepSize, surface, transform);
    }
    
    return surface;
}

void drawSurface(const Surface &surface, bool shaded)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    if (shaded)
    {
        // This will use the current material color and light
        // positions.  Just set these in drawScene();
        glEnable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // This tells openGL to *not* draw backwards-facing triangles.
        // This is more efficient, and in addition it will help you
        // make sure that your triangles are drawn in the right order.
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
    }
    else
    {        
        glDisable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        
        glColor4f(0.4f,0.4f,0.4f,1.f);
        glLineWidth(1);
    }

    glBegin(GL_TRIANGLES);
    for (unsigned i=0; i<surface.VF.size(); i++)
    {
        glNormal(surface.VN[surface.VF[i][0]]);
        glVertex(surface.VV[surface.VF[i][0]]);
        glNormal(surface.VN[surface.VF[i][1]]);
        glVertex(surface.VV[surface.VF[i][1]]);
        glNormal(surface.VN[surface.VF[i][2]]);
        glVertex(surface.VV[surface.VF[i][2]]);
    }
    glEnd();

    glPopAttrib();
}

void drawNormals(const Surface &surface, float len)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDisable(GL_LIGHTING);
    glColor4f(0,1,1,1);
    glLineWidth(1);

    glBegin(GL_LINES);
    for (unsigned i=0; i<surface.VV.size(); i++)
    {
        glVertex(surface.VV[i]);
        glVertex(surface.VV[i] + surface.VN[i] * len);
    }
    glEnd();

    glPopAttrib();
}

void outputObjFile(ostream &out, const Surface &surface)
{
    
    for (unsigned i=0; i<surface.VV.size(); i++)
        out << "v  "
            << surface.VV[i][0] << " "
            << surface.VV[i][1] << " "
            << surface.VV[i][2] << endl;

    for (unsigned i=0; i<surface.VN.size(); i++)
        out << "vn "
            << surface.VN[i][0] << " "
            << surface.VN[i][1] << " "
            << surface.VN[i][2] << endl;

    out << "vt  0 0 0" << endl;
    
    for (unsigned i=0; i<surface.VF.size(); i++)
    {
        out << "f  ";
        for (unsigned j=0; j<3; j++)
        {
            unsigned a = surface.VF[i][j]+1;
            out << a << "/" << "1" << "/" << a << " ";
        }
        out << endl;
    }
}
