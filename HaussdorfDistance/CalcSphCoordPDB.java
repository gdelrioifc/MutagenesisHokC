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

import java.io.*;
import java.util.*;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Group;

public class CalcSphCoordPDB
{
 public static LinkedList loadHydrophobs() throws Exception
  {
   LinkedList result = new LinkedList();

   result.add("ALA");
   result.add("CYS");
   result.add("PHE");
   result.add("GLY");
   result.add("ILE");
   result.add("LEU");
   result.add("MET");
   result.add("SER");
   result.add("THR");
   result.add("VAL");
   result.add("TRP");

   return result;
  }// end loadHydrophobs

 public static LinkedList loadKeysForTable(String f) throws Exception
  {
   LinkedList result = new LinkedList();

   String line="", v1="", v2="", edge="", edger="";
   StringTokenizer st = null;
   BufferedReader infile = new BufferedReader(new FileReader(f));

   while((line=infile.readLine())!=null)
    {
     st = new StringTokenizer(line,"\t");
     if(st.countTokens()==2)
      {
       v1=(String)st.nextElement();
       v2=(String)st.nextElement();

       edge=v1+" "+v2; edger=v2+" "+v1;
       if(!result.contains(edge) && !result.contains(edger)) result.add(edge);
      }
     while(st.hasMoreElements()) st.nextElement();
    }
   infile.close();

   return result;
  }// end loadKeysForTable

 public static Hashtable getTable(String f) throws Exception
  {
   Hashtable result = new Hashtable();

   int i=0;
   String token="", v1="", v2="";
   StringTokenizer st = null;
   LinkedList edges = loadKeysForTable(f);
System.err.println("No. edges="+edges.size());
   LinkedList neighbors = null;
   BufferedReader infile = new BufferedReader(new FileReader(f));

   for(i=0; i<edges.size(); i++)
    {
     token=(String)edges.get(i);
     st = new StringTokenizer(token);
     if(st.countTokens()==2)
      {
       v1=(String)st.nextElement();
       v2=(String)st.nextElement();

       if(!result.containsKey(v1))
        {
         neighbors = new LinkedList();
         neighbors.add(v2);
         result.put(v1,neighbors);
        }
       else
        {
         neighbors=(LinkedList)result.get(v1);
         if(!neighbors.contains(v2)) neighbors.add(v2);
         result.remove(v1);
         result.put(v1,neighbors);
        }
      }
     while(st.hasMoreElements()) st.nextElement();
    }

   return result;
  }// end getTable

 public static double getDistance(double x1, double x2, double y1, double y2, double z1, double z2) throws Exception
  {
   double result=-1.0;

   result=Math.pow(Math.abs(x1-x2),2.0)+Math.pow(Math.abs(y1-y2),2.0)+Math.pow(Math.abs(z1-z2),2.0);
   result=Math.sqrt(result);

   return result;
  }

 public static void printByGroup(Hashtable data4Print) throws Exception
  {
   int i=0, j=0, k=0;
   double rho=0.0, phi=0.0, theta=0.0, angulo=0.0, avgRho=0.0, avgPhi=0.0, avgTheta=0.0;
   String centerGroupRes="", centerGroupAA="", token="",rhoS1="", phiS1="", thetaS1="";
   StringTokenizer st = null;
   Hashtable center_avgValRho = new Hashtable();
   Hashtable center_avgValPhi = new Hashtable();
   Hashtable center_avgValTheta = new Hashtable();
   LinkedList rhoValues = null, phiValues=null, thetaValues=null, values=null;
   Enumeration e = null;

System.err.println("data4Print size="+data4Print.size());

   for(e = data4Print.keys(); e.hasMoreElements();)
    {
     centerGroupRes=(String)e.nextElement();
     centerGroupAA=centerGroupRes.substring(0,3);
     values=(LinkedList)data4Print.get(centerGroupRes);
     for(i=0; i<values.size(); i++)
      {
       token=(String)values.get(i);
       st = new StringTokenizer(token,"\t");
       if(st.countTokens()==3)
        {
         rhoS1=(String)st.nextElement();
         phiS1=(String)st.nextElement();
         thetaS1=(String)st.nextElement();

         if(!center_avgValRho.containsKey(centerGroupAA))
          {
//System.err.println(centerGroupAA+":"+rhoS1+":"+phiS1+":"+thetaS1);
           rhoValues = new LinkedList();
           rhoValues.add(rhoS1);
           center_avgValRho.put(centerGroupAA, rhoValues);
           phiValues = new LinkedList();
           phiValues.add(phiS1);
           center_avgValPhi.put(centerGroupAA, phiValues);
           thetaValues = new LinkedList();
           thetaValues.add(thetaS1);
           center_avgValTheta.put(centerGroupAA, thetaValues);
          }
         else
          {
//System.err.println("*"+centerGroupAA+":"+rhoS1+":"+phiS1+":"+thetaS1);
           rhoValues=(LinkedList)center_avgValRho.get(centerGroupAA);
           rhoValues.add(rhoS1);
           center_avgValRho.remove(centerGroupAA);
           center_avgValRho.put(centerGroupAA, rhoValues);
           phiValues=(LinkedList)center_avgValPhi.get(centerGroupAA);
           phiValues.add(phiS1);
           center_avgValPhi.remove(centerGroupAA);
           center_avgValPhi.put(centerGroupAA, phiValues);
           thetaValues=(LinkedList)center_avgValTheta.get(centerGroupAA);
           thetaValues.add(thetaS1);
           center_avgValTheta.remove(centerGroupAA);
           center_avgValTheta.put(centerGroupAA, thetaValues);
          }
        }
       while(st.hasMoreElements()) st.nextElement();
      }
    }

//System.err.println(center_avgValRho.size()+":"+center_avgValPhi.size()+":"+center_avgValTheta.size());
//System.err.println(center_avgValRho.toString());

   i=0;
   System.out.println("Central Residue\t<Rho>\t<Phi>\t<Theta>");
   for(e=center_avgValRho.keys(); e.hasMoreElements();)
    {
     centerGroupAA=(String)e.nextElement();
     avgRho=0.0; avgPhi=0.0; avgTheta=0.0;
     rhoValues=(LinkedList)center_avgValRho.get(centerGroupAA);
//System.err.println(centerGroupAA+":"+rhoValues.size());
      i=i+rhoValues.size();
     for(j=0; j<rhoValues.size(); j++)
      {
       rho=Double.parseDouble((String)rhoValues.get(j));
       avgRho=avgRho+rho;
      }
     avgRho=avgRho/(1.0*rhoValues.size());

     phiValues=(LinkedList)center_avgValPhi.get(centerGroupAA);
     for(j=0; j<phiValues.size(); j++)
      {
       phi=Double.parseDouble((String)phiValues.get(j));
       avgPhi=avgPhi+phi;
      }
     avgPhi=avgPhi/(1.0*phiValues.size());

     thetaValues=(LinkedList)center_avgValTheta.get(centerGroupAA);
     for(j=0; j<thetaValues.size(); j++)
      {
       theta=Double.parseDouble((String)thetaValues.get(j));
       avgTheta=avgTheta+theta;
      }
     avgTheta=avgTheta/(1.0*thetaValues.size());

     System.out.println(centerGroupAA+"\t"+avgRho+"\t"+avgPhi+"\t"+avgTheta);
    }
System.err.println("No. edges data4Print="+i);
  }// end printByGroup

 public static void main(String[] args) throws Exception
  {
   if(args.length==5)
    {
     int i=0;
     String centerRes="", neighborRes="", centerAA="", centerResNum="", neighborAA="", neighborResNum="";
     String hydrophobs=args[3].toUpperCase(), group=args[4].toUpperCase(), centerGroupRes="", centerGroupAA="", token="";
     Hashtable pdb_graph = getTable(args[1]);
     Hashtable data4Print = new Hashtable();
     LinkedList neighbors=null, values = null;
     LinkedList hydrophobsF1 = loadHydrophobs();
     Enumeration e = null;

     PDBFileReader pdbreader = new PDBFileReader();
     Structure struc = pdbreader.getStructure(args[0]);

     Chain chain = struc.getChainByPDB(args[2]);
     Group centerGroup=null, neighborGroup=null;
     AminoAcid aa1=null, aa2=null;
     Atom centerCA=null, neighborCA=null;
     double centerX=0.0, centerY=0.0, centerZ=0.0, neighborX=0.0, neighborY=0.0, neighborZ=0.0;
     double rho=0.0, phi=0.0, theta=0.0, angulo=0.0;


     if(group.equals("FALSE")) System.out.println("Central Residue\tNeighbor Residue\tRho\tPhi\tTheta");
     for(e = pdb_graph.keys(); e.hasMoreElements();)
      {
       centerRes=(String)e.nextElement();
       centerAA=centerRes.substring(0,3);
       centerResNum=centerRes.substring(3);
       centerGroup=chain.getGroupByPDB(centerResNum);

       if(centerGroup instanceof AminoAcid) 
        {
         aa1=(AminoAcid)centerGroup;
         centerCA=aa1.getCA();
         centerX=centerCA.getX();
         centerY=centerCA.getY();
         centerZ=centerCA.getZ();
         neighbors=(LinkedList)pdb_graph.get(centerRes);
       
         for(i=0; i<neighbors.size(); i++)
          {
           neighborRes=(String)neighbors.get(i);
           neighborAA=neighborRes.substring(0,3);
           neighborResNum=neighborRes.substring(3);
           neighborGroup=chain.getGroupByPDB(neighborResNum);
         
           if(neighborGroup instanceof AminoAcid)
            {
             aa2=(AminoAcid)neighborGroup;
             neighborCA=aa2.getCA();
             neighborX=neighborCA.getX();
             neighborY=neighborCA.getY();
             neighborZ=neighborCA.getZ();

             rho=getDistance(centerX, neighborX, centerY, neighborY, centerZ, neighborZ);

             phi=Math.acos((neighborZ-centerZ)/rho);
//System.err.print(rho+"\t"+phi+"\t");
             angulo=((phi)/Math.PI)*180.0;
//System.err.print(rho+"\t"+angulo+"\t");
             phi=angulo;

             theta=Math.atan(Math.abs(neighborY-centerY)/Math.abs(neighborX-centerX));
//System.err.println(theta);
             angulo=(theta/Math.PI)*180.0;
//System.err.println(angulo);
//             if(angulo<0) angulo=180+angulo;
             theta=angulo;

             if(hydrophobs.equals("TRUE") & hydrophobsF1.contains(centerAA) && hydrophobsF1.contains(neighborAA))
              {
               if(group.equals("TRUE")) 
                {
                 if(!data4Print.containsKey(centerRes))
                  {
                   values = new LinkedList();
                   values.add(String.valueOf(rho+"\t"+phi+"\t"+theta));
                   data4Print.put(centerRes,values);
                  }
                 else
                  {
                   values=(LinkedList)data4Print.get(centerRes);
                   values.add(String.valueOf(rho+"\t"+phi+"\t"+theta));
                   data4Print.remove(centerRes);
                   data4Print.put(centerRes, values);
                  }
                }
               else System.out.println(centerRes+"\t"+neighborRes+"\t"+rho+"\t"+phi+"\t"+theta);
              }
             else if(hydrophobs.equals("FALSE"))
              {
               if(group.equals("TRUE"))
                {
                 if(!data4Print.containsKey(centerRes))
                  {
                   values = new LinkedList();
                   values.add(String.valueOf(rho+"\t"+phi+"\t"+theta));
                   data4Print.put(centerRes,values);
                  }
                 else
                  {
                   values=(LinkedList)data4Print.get(centerRes);
                   values.add(String.valueOf(rho+"\t"+phi+"\t"+theta));
                   data4Print.remove(centerRes);
                   data4Print.put(centerRes, values);
                  }
                }
               else System.out.println(centerRes+"\t"+neighborRes+"\t"+rho+"\t"+phi+"\t"+theta);
              }
            }
          }
        }
      }

     if(group.equals("TRUE")) printByGroup(data4Print);
    }
   else
    {
     System.err.println("Gnomovision version 69, Copyright (C) 2016 Gabriel Del Rio");
     System.err.println("Gnomovision comes with ABSOLUTELY NO WARRANTY."); 
     System.err.println("Please read the conditions at http://www.gnu.org/licenses/gpl-2.0.html.\n\n");

     System.err.println("Usage:\njava CalcSphCoordPDB <pdb_file> <pdb_graph> <chain> <hydrophobs> <group>");
     System.err.println("<pdb_file> a file in PDB format; e.g. 1hiv.pdb in this dir.");
     System.err.println("<pdb_graph> a 2-column tab-delimited file generated by program doGraphFromPDB in this dir; e.g., 1hiv.graph");
     System.err.println("<chain> the name of the chain in the <pdb_file> that was used to obtain <pdb_graph>; e.g., A");
     System.err.println("<hydrophobs> is true or false.");
     System.err.println("<group> is true or false.");
     System.err.println("The program will compute the spherical coordinates (rho, phi and theta) for each pair of residues");
     System.err.println("listed in <pdb_graph> derived from the cartesians coordinates in <pdb_file>.");
     System.err.println("If <hydrophobs>==true then, only hydrophobic residues will be considered in this analisis.");
     System.err.println("If <group>==true then, 20 groups, each per every aa, will be created and the average angles and distances will be reported.");
     System.err.println("The output will be sent to the standard output.");
     System.err.println("");
     System.err.println("");
    }
  }
}
