import java.io.*;
import java.util.*;

public class GetCorrelatedMutationsAtClosePositions
{
 public static Hashtable tableAA;
 public static LinkedList protsID;

 public static void loadTableAA() throws Exception
  {
   tableAA = new Hashtable();

   tableAA.put("ALA","A");
   tableAA.put("CYS","C");
   tableAA.put("ASP","D");
   tableAA.put("GLU","E");
   tableAA.put("PHE","F");
   tableAA.put("GLY","G");
   tableAA.put("HIS","H");
   tableAA.put("ILE","I");
   tableAA.put("LYS","K");
   tableAA.put("LEU","L");
   tableAA.put("MET","M");
   tableAA.put("ASN","N");
   tableAA.put("PRO","P");
   tableAA.put("GLN","Q");
   tableAA.put("ARG","R");
   tableAA.put("SER","S");
   tableAA.put("THR","T");
   tableAA.put("VAL","V");
   tableAA.put("TRP","W");
   tableAA.put("TYR","Y");

  }// end loadTableAA

 public static Hashtable loadHSSP(String f) throws Exception
  {
   Hashtable result = new Hashtable();

   protsID = new LinkedList();
   String line="", aa="", pos="", res="", token="";
   BufferedReader infile = new BufferedReader(new FileReader(f));

   breakable:
   while((line=infile.readLine())!=null)
    {
     if(line.startsWith("## PROTEINS : identifier and alignment statistics"))
      {
       line=infile.readLine();
       while(!(line=infile.readLine()).startsWith("## ALIGNMENTS"))
       {
        token=line.substring(8,20).trim();
//System.err.println(token);
        protsID.add(token);
       }
      }
     else if(line.indexOf(" SeqNo  PDBNo AA")==0)
      {
       while((line=infile.readLine()).indexOf("## ALIGNMENTS")==-1)
        {
//System.err.println(line);
         if(line.startsWith("## SEQUENCE PROFILE AND ENTROPY"))
          {
           break breakable;
          }
         else
          {
           pos=line.substring(7,11).trim();
           aa=line.substring(13,15).trim().toUpperCase();
//System.err.println(aa+":"+pos);

//if(aa.equals("C") && pos.equals("69")) System.err.println(aa+pos+":"+token);
           if(aa.matches("[A-Y]") && !pos.equals(""))
            {
             res=aa+pos;
             token=line.substring(51).toUpperCase();
//System.err.println(line);
//System.err.println("Mira="+res+":"+token);

             if(!result.containsKey(res)) result.put(res,token);
             else
              {
               token=token+(String)result.get(res);
               result.remove(res);
               result.put(res,token);
              }
            }
          }
        }
      }
    }
   infile.close();

   return result;
  }// end loadHSSP

 public static Hashtable loadTable(String f) throws Exception
  {
   Hashtable result = new Hashtable();

   String line="", token="", r1="", r2="", aa="", pos="";
   StringTokenizer st1 = null, st2=null;
   LinkedList contacts = null;
   BufferedReader infile = new BufferedReader(new FileReader(f));

   while((line=infile.readLine())!=null)
    {
     st1 = new StringTokenizer(line,"\t");
     if(st1.countTokens()==2)
      {
       r1=((String)st1.nextElement()).trim();
       aa=(String)tableAA.get(r1.substring(0,3));
       pos=r1.substring(3);
       r1=aa+pos;
       token=(String)st1.nextElement();

       st2 = new StringTokenizer(token,",");
       contacts = new LinkedList();
       while(st2.hasMoreElements())
        {
         r2=((String)st2.nextElement()).trim();
         aa=(String)tableAA.get(r2.substring(0,3));
         pos=r2.substring(3);
         r2=aa+pos;
         contacts.add(r2);
        }
      }
     while(st1.hasMoreElements()) st1.nextElement();

//System.err.println(r1+":"+contacts.toString());
     result.put(r1,contacts);
    }
   infile.close();

   return result;
  }// end loadTable

 public static LinkedList getVariantPositions(String aa1, String variants) throws Exception
  {
   LinkedList result = new LinkedList();

   int i=0;
   String r2="", aa2="";

//System.err.println("* "+aa1+":"+variants);
   if(variants!=null)
    {

     for(i=0; i<variants.length(); i++)
      {
       r2=String.valueOf(variants.charAt(i));
       aa2=r2.substring(0,1);
       if(aa1.equals(aa2)) result.add("0");
       else result.add("1");
      }
    }

   return result;
  }// end getVariantPositions

 public static LinkedList getCorrelatedMutations(LinkedList pos_var1, LinkedList pos_var2) throws Exception
  {
   LinkedList result = new LinkedList();

   int i=0, j=0;
   String v1="", v2="";

   for(i=0; i<pos_var1.size(); i++)
    {
     v1=(String)pos_var1.get(i);
     v2=(String)pos_var2.get(i);

     if(v1.equals("1") && v2.equals("1")) result.add(i,"1");
     else result.add(i,"0");
    }

   return result;
  }// end getCorrelatedMutations

 public static String printCorMutData(LinkedList cor_muts, String var1, String var2) throws Exception
  {
   String result="No correlation found";

   int i=0, index=0, count=0;
   String pid="";

   for(i=0; i<cor_muts.size(); i++)
    {
     index=Integer.parseInt((String)cor_muts.get(i));
     if(index==1)
      {
       count=count+1;
       if(result.equals("No correlation found"))
        {
         pid=(String)protsID.get(i);
         result=pid+"("+String.valueOf(var1.charAt(i))+":"+String.valueOf(var2.charAt(i))+")";
        }
       else
        {
         pid=(String)protsID.get(i);
         result=result+","+pid+"("+String.valueOf(var1.charAt(i))+":"+String.valueOf(var2.charAt(i))+")";
        }
      }
    }

   if(var1.length()!=var2.length()) result = result + "\t" + String.valueOf(count)+" / "+var1.length()+" OR "+var2.length();
   else result = result + "\t" + String.valueOf(count)+" / "+var1.length();

   return result;
  }// end printCorMutData

 public static void main(String[] args) throws Exception
  {
   if(args.length==2)
    {
     int i=0;
     String r1="", aa1="", r2="", aa2="", var1="", var2="";
     LinkedList contacts = null, pos_var1=null, pos_var2=null, cor_muts=null, wo_variants=new LinkedList();
     StringTokenizer st = null;
     loadTableAA();
     Hashtable residue_variance = loadHSSP(args[0]);
     Hashtable residue_contacts = loadTable(args[1]);

     System.out.println("Correlated positions\tProt IDs\tCount of correlations");
     for(Enumeration e = residue_contacts.keys(); e.hasMoreElements();)
      {
       r1=(String)e.nextElement();
       aa1=r1.substring(0,1);
       contacts=(LinkedList)residue_contacts.get(r1);
       var1=(String)residue_variance.get(r1);
//System.err.println("*:"+r1+":"+var1);
       if(var1!=null)
        {
         pos_var1=getVariantPositions(aa1,var1);

         for(i=0; i<contacts.size(); i++)
          {
           r2=(String)contacts.get(i);
           aa2=r2.substring(0,1);
           var2=(String)residue_variance.get(r2);
//System.err.println(r1+":"+r2+":"+var2);
           if(var2!=null)
            {
             pos_var2=getVariantPositions(aa2,var2);

             if(pos_var1.size()==pos_var2.size()) 
              {
               cor_muts=getCorrelatedMutations(pos_var1, pos_var2);
               System.out.println(r1+":"+r2+"\t"+printCorMutData(cor_muts, var1, var2));
              }
             else
              {
               System.out.println("* "+r1+":"+r2+"\t"+"Number of variant positions is different\t0 / "+pos_var1.size()+" OR "+pos_var2.size());
              }
            }
           else if(!wo_variants.contains(r2)) wo_variants.add(r2);
          }
         System.out.println("");
        }
       else if(!wo_variants.contains(r1)) wo_variants.add(r1);
      }
     System.out.println("Residues in "+args[1]+" without variants in "+args[0]+":\t"+wo_variants.toString());
    }
   else
    {
     System.err.println("Usage:\njava GetCorrelatedMutationsAtClosePositions <hssp_file> <pdb_contacts>");
     System.err.println("<hssp_file> name of an HSSP file; e.g., 4tsv.hssp");
     System.err.println("<pdb_contacts> file generated by program GetContacts; e.g. 4tsv.graph_contacts");
     System.err.println("This program will print out the residues in <hssp_file> that are in contact in <pdb_contacts> for the reference protein");
     System.err.println("(e.g., 4tsv). Then, it will identify simultaneous substitutions at these positions in any of the protein");
     System.err.println("sequences aligned in <hssp_file>. Those positions that are in contacto and simultaneously changed are considered");
     System.err.println("compensatory/correlated. The frequency of these correlations is reported too.");
     System.err.println("The output will be sent to the standard output.");
     System.err.println("");
     System.err.println("");
    }
  }
}
