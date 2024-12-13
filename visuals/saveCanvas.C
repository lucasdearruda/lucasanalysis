void sC(TCanvas *myCanvas, string path_to_save = Form("%s/images/",getenv("PWD"))) {
    // Obtém a data atual e o horário
    time_t now = time(0);
    tm *ltm = localtime(&now);
    
    int year = 1900 + ltm->tm_year;
    int month = 1 + ltm->tm_mon;
    int day = ltm->tm_mday;
    int hour = ltm->tm_hour;
    int minute = ltm->tm_min;
    
    // Formata a data
    string date = Form("%04d-%02d-%02d", year, month, day);
    
    // Verifica o próximo número disponível
    int nextNumber = 1;
    void *dir = gSystem->OpenDirectory(path_to_save.c_str());
    const char *file;
    while ((file = gSystem->GetDirEntry(dir))) {
        if (string(file).find(date) != string::npos) {
            int number = atoi(string(file).substr(11, 2).c_str());
            if (number >= nextNumber) {
                nextNumber = number + 1;
            }
        }
    }
    gSystem->FreeDirectory(dir);
    
    // Formata o nome do arquivo
    string filename = Form("%s%s-%02d.root", path_to_save.c_str(), date.c_str(), nextNumber);
    
    // Escreve o nome do arquivo e a hora no canto inferior esquerdo do canvas
    myCanvas->cd();
    
    TPad *newpad=new TPad("newpad","a transparent pad",0,0,1,1);
    newpad->SetFillStyle(4000);
    newpad->Draw();
    newpad->cd();
    //test:
    myCanvas->SetEditable(kTRUE);
    TLatex text;
    text.SetTextSize(0.03);
    text.SetTextColor(kRed);
    text.DrawLatexNDC(0.05, 0.93, Form("%s,  %02d:%02d", filename.c_str(), hour, minute));
    

    // Salva o canvas
    myCanvas->SaveAs(filename.c_str());
    filename = Form("%s%s-%02d.png", path_to_save.c_str(), date.c_str(), nextNumber);
    myCanvas->SaveAs(filename.c_str());
    
}

