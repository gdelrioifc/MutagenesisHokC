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

public class AnalyzeCorrMuts
{
 public static Hashtable getRatios(String f) throws Exception
  {
   Hashtable result = new Hashtable();

   double n1=0.0, n2=0.0, ratio=0.0;
   String line="", par="";
   String[] tokens = null;
   BufferedReader infile = new BufferedReader(new FileReader(f));

   while((line=infile.readLine())!=null)
    {
     if(!line.startsWith("Correlated positions") && !line.startsWith("Residues in") && !line.equals(""))
      {
       tokens = line.split("[ \t]");
       par=tokens[0];
//System.err.println(line);
//System.err.println(tokens[tokens.length-3]);
//System.err.println(tokens[tokens.length-1]);

       n1=Double.parseDouble(tokens[tokens.length-3]);
       n2=Double.parseDouble(tokens[tokens.length-1]);
       ratio=n1/n2;

       if(!result.containsKey(par)) result.put(par,String.valueOf(ratio));
       else
        {
         ratio=ratio+Double.parseDouble((String)result.get(par));
         result.remove(par);
         result.put(par, String.valueOf(ratio));
        }
      }
    }
   infile.close();

   return result;
  }// end getRatios

 public static Hashtable getSustitutionCounts(String f) throws Exception
  {
   Hashtable result = new Hashtable();

   int i=0, count=0;
   String line="", par="", par1="", par2="", s1="", s2="", token1="", token2="";
   String[] tokens1 = null, tokens2 = null;
   BufferedReader infile = new BufferedReader(new FileReader(f));

   while((line=infile.readLine())!=null)
    {
     if(!line.startsWith("Correlated positions") && !line.startsWith("Residues in") && !line.equals(""))
      {
       tokens1 = line.split("\t");
       par=tokens1[0];
       par1=par.substring(0,1);
       par2=par.substring(par.indexOf(":")+1, par.indexOf(":")+2);
//System.err.println(par+":"+par1+"-"+par2);
       token1=tokens1[1];

       if(!token1.equals("No correlation found"))
        {
         tokens2 = token1.split(",");

         for(i=0; i<tokens2.length; i++)
          {
//System.err.println("*"+line);
           s1=tokens2[i].substring(tokens2[i].indexOf("(")+1, tokens2[i].indexOf(":"));
           s2=tokens2[i].substring(tokens2[i].indexOf(":")+1, tokens2[i].indexOf(")"));
//System.err.println(tokens2[i]+":"+s1+"-"+s2);

           token2=par1+"-"+s1+":"+par2+"-"+s2;
           if(!result.containsKey(token2)) result.put(token2,"1");
           else
            {
             count=Integer.parseInt((String)result.get(token2)) + 1;
             result.remove(token2);
             result.put(token2,String.valueOf(count));
            }
          }
        }
      }
    }
   infile.close();

   return result;
  }// end getSustitutionCounts

 public static void main(String[] args) throws Exception
  {
   if(args.length==2)
    {
     int i=0, count=0;
     String token="", value="", flag=args[1];
     Enumeration e = null;

     if(flag.equalsIgnoreCase("-b") || flag.equalsIgnoreCase("-fp"))
      {
       Hashtable par_ratio = getRatios(args[0]);

       System.out.println("Index\tAA in contact\tFrequency");
       for(e=par_ratio.keys(); e.hasMoreElements();)
        {
         token=(String)e.nextElement();
         value=(String)par_ratio.get(token);
         i=i+1;
         if(!value.equals("0.0")) count=count+1;
         System.out.println(i+"\t"+token+"\t"+value);
        }
       System.out.println("No. of different paired aas found="+par_ratio.size());
       System.out.println("No. of paired aas with correlated mutations="+count);
       System.out.println("****");
      }

     if(flag.equalsIgnoreCase("-b") || flag.equalsIgnoreCase("-sp"))
      {
       Hashtable sustitution_count = getSustitutionCounts(args[0]);

       i=0; count=0;
       System.out.println("Index\tCorrelated mutation in contact\tCount");
       for(e=sustitution_count.keys(); e.hasMoreElements();)
        {
         token=(String)e.nextElement();
         value=(String)sustitution_count.get(token);
         i=i+1;
         count=count+Integer.parseInt(value);
         System.out.println(i+"\t"+token+"\t"+value);
        }
       System.out.println("No. of different correlated mutations found="+sustitution_count.size());
       System.out.println("No. of correlated mutations found="+count);
      }
    }
   else
    {
     System.err.println("Gnomovision version 69, Copyright (C) 2016 Gabriel Del Rio");
     System.err.println("Gnomovision comes with ABSOLUTELY NO WARRANTY."); 
     System.err.println("Please read the conditions at http://www.gnu.org/licenses/gpl-2.0.html.\n\n");

     System.err.println("Uso:\njava AnalyzeCorrMuts <infile> <flag>");
     System.err.println("<infile> archivo de salida del programa ../HSSP/GetCorrelatedMutationsAtClosePositions; e.g., 4tsv.getCorrMuts");
     System.err.println("<flag> -fp, -sp or -b, see below.");
     System.err.println("El programa reportara la frecuencia normalizada de cada par de aas en");
     System.err.println("contacto (<flag>=-fp), la del numero de sustituciones correlacionadas encontradas por");
     System.err.println("par de aa en contacto (<flag>=-sp) o ambas (<flag>=-b).");
     System.err.println("La salida se hara al standard output.");
     System.err.println("");
     System.err.println("");
    }
  }
}
