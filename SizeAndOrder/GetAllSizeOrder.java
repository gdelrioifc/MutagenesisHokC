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

public class GetAllSizeOrder
{
 public static String getData(String f) throws Exception
  {
   String result=f.substring(f.lastIndexOf(File.separator)+1,f.indexOf("."));

   String line="", size="", order="";
   BufferedReader infile = new BufferedReader(new FileReader(f));

   while((line=infile.readLine())!=null)
    {
     if(line.startsWith("Size=")) size=line.substring(6);
     if(line.startsWith("Order=")) order=line.substring(7);
    }
   infile.close();

   result=result+" "+size+" "+order;
   return result;
  }// end getData

 public static void main(String[] args) throws Exception
  {
   if(args.length==1)
    {
     String line="", pdbChain="", size="", order="", token="";
     StringTokenizer st = null;
     BufferedReader infile = new BufferedReader(new FileReader(args[0]));

     while((line=infile.readLine())!=null)
      {
       token=getData(line.trim());
       st = new StringTokenizer(token);
       if(st.countTokens()==3)
        {
         pdbChain=(String)st.nextElement();
         size=(String)st.nextElement();
         order=(String)st.nextElement();

         System.out.println(pdbChain+"\t"+size+"\t"+order);
        }
       while(st.hasMoreElements()) st.nextElement();
      }
     infile.close();
    }
   else
    {
     System.err.println("Usage:\njava GetAllSizeOrder <infile>");
     System.err.println("<infile> a file with a list of file names obtained with the program calcSizeOrderGraph.");
     System.err.println("The program will print to the standard out the pdb-chain, size and order reported for all the files in <infiel>.");
     System.err.println("");
     System.err.println("");
    }
  }
}
