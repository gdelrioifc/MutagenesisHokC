/*
Copyright (C) 2016 Gabriel Del Rio

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/


import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.algorithm.distance.DiscreteHausdorffDistance;
import com.vividsolutions.jts.geom.util.AffineTransformation;
import com.vividsolutions.jts.algorithm.Centroid;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import java.lang.Double;

class GeometryHandler{
	public static void main(String args[]){
		GeometryHandler f = new GeometryHandler();
		Geometry g1 = f.toGeometry(f.readCoordinate(args[0]));
		Geometry g2 = f.toGeometry(f.readCoordinate(args[1]));
		double hdrfDist;
		hdrfDist = DiscreteHausdorffDistance.distance(g1,g2);
		System.out.println(hdrfDist);
		g1 = f.rotate(g1, 30);
		hdrfDist = DiscreteHausdorffDistance.distance(g1,g2);
		System.out.println(hdrfDist);
		g2 = f.translate(g2, 2,4);
		hdrfDist = DiscreteHausdorffDistance.distance(g1,g2);
		System.out.println(hdrfDist);
	}

	public Coordinate[] readCoordinate(String filename){
		ArrayList<Coordinate> cl = new ArrayList<Coordinate>();
    	try{
			BufferedReader br = new BufferedReader(new FileReader(filename));
			br.readLine();
			String line;
			String[] colum=null;
			while ((line = br.readLine()) != null) {
				colum = line.split("\t");
				//System.out.println(colum[0]+" "+ colum[1]+" "+colum[2]+" "+colum[3]+" "+colum[4]);
				Coordinate c = new Coordinate();
				c.x=Double.parseDouble(colum[3]);
				c.y=Double.parseDouble(colum[4]);
				//c.z=Double.parseDouble(colum[4]);
				cl.add(c);
    		}
    	} catch (Exception e) {
    		e.printStackTrace();
    	}
		Coordinate[] c = new Coordinate[cl.size()]; 
		cl.toArray(c);
		return c;
	}

	public Geometry toGeometry(Coordinate[] c){
		CoordinateArraySequence cas = new CoordinateArraySequence(c);
		Envelope e = new Envelope();
		e = cas.expandEnvelope(e);
		GeometryFactory gf = new GeometryFactory();
		Geometry g = gf.toGeometry(e);
		return g;
	}

	public Geometry rotate(Geometry g, double theta){
		AffineTransformation trans = new AffineTransformation();
		Centroid cent = new Centroid(g);
		Coordinate c = cent.getCentroid();
		trans = trans.rotate(theta, c.x, c.y);
		Geometry gOut = trans.transform(g);
		return gOut;
	}

	public Geometry translate(Geometry g, double moveX, double moveY){
		AffineTransformation trans = new AffineTransformation();
		trans = trans.translate(moveX, moveY);
		Geometry gOut = trans.transform(g);
		return gOut;
	}
}
