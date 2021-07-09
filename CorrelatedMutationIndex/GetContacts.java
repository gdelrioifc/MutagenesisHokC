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

public class GetContacts
{
 public static void main(String[] args) throws Exception
  {
   if(args.length==1)
    {
     String line="", r1="", r2="";
     StringTokenizer st = null;
     Hashtable residue_contacts = new Hashtable();
     LinkedList contacts = null;
     BufferedReader infile = new BufferedReader(new FileReader(args[0]));

     while((line=infile.readLine())!=null)
      {
       st = new StringTokenizer(line,"\t");
       if(st.countTokens()==2)
        {
         r1=(String)st.nextElement();
         r2=(String)st.nextElement();

         if(!residue_contacts.containsKey(r1))
          {
           contacts = new LinkedList();
           contacts.add(r2);
           residue_contacts.put(r1,contacts);
          }
         else
          {
           contacts=(LinkedList)residue_contacts.get(r1);
           if(!contacts.contains(r2)) contacts.add(r2);
           residue_contacts.remove(r1);
           residue_contacts.put(r1,contacts);
          }

/*         
         if(!residue_contacts.containsKey(r2))
          {
           contacts = new LinkedList();
           contacts.add(r1);
           residue_contacts.put(r2,contacts);
          }
         else
          {
           contacts=(LinkedList)residue_contacts.get(r2);
           if(!contacts.contains(r1)) contacts.add(r1);
           residue_contacts.remove(r2);
           residue_contacts.put(r2,contacts);
          }
*/
        }
       while(st.hasMoreElements()) st.nextElement();
      }
     infile.close();

     for(Enumeration e=residue_contacts.keys(); e.hasMoreElements();)
      {
       r1=(String)e.nextElement();
       contacts=(LinkedList)residue_contacts.get(r1);
       r2=(contacts.toString()).replaceAll("\\[","");
       r2=r2.replaceAll("\\]","");
       System.out.println(r1+"\t"+r2);
      }
    }
   else
    {
     System.err.println("Gnomovision version 69, Copyright (C) 2016 Gabriel Del Rio");
     System.err.println("Gnomovision comes with ABSOLUTELY NO WARRANTY."); 
     System.err.println("Please read the conditions at http://www.gnu.org/licenses/gpl-2.0.html.\n\n");

     System.err.println("Uso:\njava GetContacts <infile>");
     System.err.println("<infile> archivo generado por doGraphFromPDB; e.g., 4tsv.graph");
     System.err.println("El programa transformara el <infile> en un formato leible por");
     System.err.println("el programa GetCorrelatedMutationsAtClosePositions.");
     System.err.println("El output se hara al standard output.");
     System.err.println("");
     System.err.println("");
    }
  }
}
