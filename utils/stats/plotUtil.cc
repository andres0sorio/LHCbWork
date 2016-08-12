#include <plotUtil.h>

int findMaximum( int ref, const std::vector<double> & vec) {

  int ans = 0;
  int size = vec.size();
  
  ans = ref;

  for(int i = 0; i < size ; i++) {
    if(vec[i] < ans) {
      //don't do anything
    } else if(vec[i] > ans) {
      ans = int (vec[i]);
    }
  }

  return ans;

}

int findMaximum( int ref, const std::vector<int> & vec) {

  int ans = 0;
  int size = vec.size();
  
  ans = ref;

  for(int i = 0; i < size ; i++) {
    if(vec[i] < ans) {
      //don't do anything
    } else if(vec[i] > ans) {
      ans = int (vec[i]);
    }
  }

  return ans;

}

void setGraphOptions(TGraph *graph, char *option1, char *option2) {
  
  graph->SetMarkerStyle(21);
  graph->SetMarkerSize(0.6);
  graph->SetMinimum(0);
  graph->SetTitle("");
  TAxis *ax1 = graph->GetXaxis();
  TAxis *ax2 = graph->GetYaxis();
  setAxeOptions(ax1);
  setAxeOptions(ax2);
  ax1->SetTitle(option1);
  ax2->SetTitle(option2);
}

void setFrameOptions(TPad *pad, char *option1, char *option2) {
  
  pad->GetFrame()->SetFillColor(10);
  TAxis *ax1 = ((TH1F*)pad->FindObject("hframe"))->GetXaxis();
  TAxis *ax2 = ((TH1F*)pad->FindObject("hframe"))->GetYaxis();

  setAxeOptions(ax1);
  setAxeOptions(ax2);

  //ax1->SetTitle(option1);
  ax1->SetNdivisions(0);
  ax2->SetTitle(option2);


}

void setAxeOptions(TAxis * ax) {
    
  ax->SetLabelSize(0.03);
  ax->SetLabelFont(42);
  ax->SetLabelOffset(0.006);
  
  ax->SetTitleSize(0.03);
  ax->SetTitleFont(42);
  ax->SetTitleOffset(1.5);
  
}


void createLabels ( TPad * pad, 
		    TPaveLabel *labels[], 
		    std::string *names,
		    int max)
{

  double x1(0.0), x2(0.0), y1(0.0), y2(0.0);
  double dx(0.0);
  double dy(0.0);
  double gap(0.0);
  int event(0);
  char buf [15];
  int side_limit = max;
  
  TLine *lines[3*max];
  
  dx = 0.08;
  dy = 0.06;
  gap = 0.0175;

  for( int i = 0; i < side_limit; i++) {
    
    y1 = (0.)*dy + (0.+0.001)*gap;
    y2 = y1+dy;
    
    x1 = (i)*dx + (i+1.)*gap;
    x2 = x1+dx;
    
    //sprintf(buf, "%d", (event+0));
    sprintf(buf, "%s", names[i].c_str());

    labels[event] = new TPaveLabel(x1,y1,x2,y2,buf,"NDC");
    setLabelOptions(labels[event]);
    labels[event]->Draw();
    joinLabels(pad , labels[event], lines, i );
    event++;
    
  }
  
  
}

void joinLabels ( TPad *pad, TPaveLabel * box, TLine **lines, int pos )
{
  
  double fixed_x(0.0), fixed_y(0.0);
  
  double x1(0.0), x2(0.0); 
  double x3(0.0), y1(0.0); 
  double y2(0.0), y3(0.0);
  double x4(0.0), y4(0.0);
  
  double offset(0.0);
  double gap(0.0);
  
  offset = 0.1;
  gap = 0.07;
  
  int i1 = pos;
  int i2 = pos+1;
  int i3 = pos+2;

  x1 = (box->GetX2()) - (0.08 / 2.0);
  y1 = box->GetY2();
  fixed_x = 0.1726 + (pos*0.07277777);
  fixed_y = 0.07;
  
  x2 = fixed_x;
  y2 = fixed_y;

  x3 = x2;
  y3 = 0.14;
  
  x4 = x2;
  y4 = 0.5;
  
  lines[i1] = new TLine();
  lines[i2] = new TLine();
  lines[i3] = new TLine();

  setLineOptions(lines[i1]);
  setLineOptions(lines[i2]);
  setLineOptions(lines[i3]);

  //lines[i3]->SetLineStyle(4);
  
  lines[i1]->DrawLineNDC(x1,y1,fixed_x,fixed_y);
  lines[i2]->DrawLineNDC(x2,y2,x3,y3);
  lines[i3]->DrawLineNDC(x3,y4,x4,y4);
}

void setLabelOptions ( TPaveLabel *label)
{
  
  label->SetFillColor(10);
  label->SetBorderSize(0);
  label->SetTextSize(0.5);
  label->SetTextAlign(22);
  //label->SetTextAngle(15.0);
  
}

void setLineOptions ( TLine *line) 
{
  
  line->SetLineStyle(1);
  line->SetLineColor(1);
  line->SetLineWidth(1);
  
}

void setCutNames(std::string *names) {

  names[0] = "M_{recoil}";
  names[1] = "P_{trans}";
  names[2] = "E_{trans}";
  names[3] = "YCutEt";
  names[4] = "Cos pm";
  names[5] = "Cos em";
  names[6] = "n_Ch.trks";
  names[7] = "E_{around}";
  names[8] = "Prob(#chi^{2})";
  names[9] = "All";
  
}
