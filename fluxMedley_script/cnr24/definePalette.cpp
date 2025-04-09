


void SetWhichPalette(int id=1) {
    Int_t nColors;  // Define number of colors
    Int_t* colors;  // Define pointer to array of colors
    
    if (id == 1) {  // Bang Wong's
        cout << "Bang Wong's Palette!" << endl;
        nColors = 9;
        Int_t bangWongColors[] = {
            TColor::GetColor("#000000"),  // Black
            TColor::GetColor("#E69F00"),  // Yellow
            TColor::GetColor("#56B4E9"),  // Cyan
            TColor::GetColor("#009E73"),  // Green
            TColor::GetColor("#F0E442"),  // Yellow
            TColor::GetColor("#0072B2"),  // Blue
            TColor::GetColor("#D55E00"),  // Orange
            TColor::GetColor("#CC79A7")   // Rose
        };
        colors = bangWongColors;
    } else if (id == 2) {  // IBM's
        cout << "IBM's Palette!" << endl;
        nColors = 5;
        Int_t ibmColors[] = {
            TColor::GetColor("#648FFF"),  // Blue
            TColor::GetColor("#785EF0"),  // Purple
            TColor::GetColor("#DC267F"),  // Salmon
            TColor::GetColor("#FE6100"),  // Orange
            TColor::GetColor("#FFB000")   // Gold
        };
        colors = ibmColors;
    } else if (id == 3) {  // Paul Tol's
        cout << "Paul Tol's Palette!" << endl;
        nColors = 8;
        Int_t paulTolColors[] = {
            TColor::GetColor("#332288"),  // Purple
            TColor::GetColor("#117733"),  // Dark Green
            TColor::GetColor("#44AA99"),  // Aquamarine
            TColor::GetColor("#88CCEE"),  // Sky Blue
            TColor::GetColor("#DDCC77"),  // LightWood
            TColor::GetColor("#CC6677"),  // Salmon
            TColor::GetColor("#AA4499"),  // Violet
            TColor::GetColor("#882255")   // Satin
        };
        colors = paulTolColors;
    } else if (id == 4) {  // Sunset Discrete
        cout << "Sunset Discrete Palette!" << endl;
        nColors = 11;
        Int_t sunsetColors[] = {
            TColor::GetColor("#364B9A"),  // Dark Blue
            TColor::GetColor("#4A7BB7"),  // Light Blue
            TColor::GetColor("#6EA6CD"),  // Sky Blue
            TColor::GetColor("#98CAE1"),  // Pale Blue
            TColor::GetColor("#C2E4EF"),  // Very Pale Blue
            TColor::GetColor("#EAECCC"),  // Cream
            TColor::GetColor("#FEDA8B"),  // Light Yellow
            TColor::GetColor("#FDB366"),  // Yellow-Orange
            TColor::GetColor("#F67E4B"),  // Orange
            TColor::GetColor("#DD3D2D"),  // Red
            TColor::GetColor("#A50026")   // Dark Red
        };
        colors = sunsetColors;
    } else if (id == 5) {  // Nightfall Discrete
        cout << "Nightfall Discrete Palette!" << endl;
        nColors = 9;
        Int_t nightfallColors[] = {
            TColor::GetColor("#125A56"),  // Dark Teal
            TColor::GetColor("#238F9D"),  // Teal
            TColor::GetColor("#60BCE9"),  // Light Blue
            TColor::GetColor("#C6DBED"),  // Very Light Blue
            TColor::GetColor("#ECEADA"),  // Light Beige
            TColor::GetColor("#F9D576"),  // Light Yellow
            TColor::GetColor("#FD9A44"),  // Orange
            TColor::GetColor("#E94C1F"),  // Red
            TColor::GetColor("#A01813")   // Dark Red
        };
        colors = nightfallColors;
    } else if (id == 6) {  // BuRd Discrete
        cout << "BuRd Discrete Palette!" << endl;
        nColors = 9;
        Int_t buRdColors[] = {
            TColor::GetColor("#2166AC"),  // Dark Blue
            TColor::GetColor("#4393C3"),  // Light Blue
            TColor::GetColor("#92C5DE"),  // Pale Blue
            TColor::GetColor("#D1E5F0"),  // Very Pale Blue
            TColor::GetColor("#F7F7F7"),  // Light Gray
            TColor::GetColor("#FDDBC7"),  // Light Peach
            TColor::GetColor("#F4A582"),  // Peach
            TColor::GetColor("#D6604D"),  // Coral
            TColor::GetColor("#B2182B")   // Dark Red
        };
        colors = buRdColors;
    } else if (id == 7) {  // PRGn Discrete
        cout << "PRGn Discrete Palette!" << endl;
        nColors = 9;
        Int_t prGnColors[] = {
            TColor::GetColor("#762A83"),  // Dark Purple
            TColor::GetColor("#9970AB"),  // Purple
            TColor::GetColor("#C2A5CF"),  // Light Purple
            TColor::GetColor("#E7D4E8"),  // Very Light Purple
            TColor::GetColor("#F7F7F7"),  // Light Gray
            TColor::GetColor("#D9F0D3"),  // Pale Green
            TColor::GetColor("#ACD39E"),  // Light Green
            TColor::GetColor("#5AAE61"),  // Green
            TColor::GetColor("#1B7837")   // Dark Green
        };
        colors = prGnColors;
    } else if (id == 8) {  // YlOrBr Discrete
        cout << "YlOrBr Discrete Palette!" << endl;
        nColors = 9;
        Int_t ylOrBrColors[] = {
            TColor::GetColor("#FFFFE5"),  // Very Pale Yellow
            TColor::GetColor("#FFF7BC"),  // Pale Yellow
            TColor::GetColor("#FEE391"),  // Light Yellow
            TColor::GetColor("#FEC44F"),  // Yellow
            TColor::GetColor("#FB9A29"),  // Orange
            TColor::GetColor("#EC7014"),  // Dark Orange
            TColor::GetColor("#CC4C02"),  // Brown
            TColor::GetColor("#993404"),  // Dark Brown
            TColor::GetColor("#662506")   // Very Dark Brown
        };
        colors = ylOrBrColors;
    } else if (id == 9) {  // WhOrBr
        cout << "WhOrBr Palette!" << endl;
        nColors = 9;
        Int_t whOrBrColors[] = {
            TColor::GetColor("#FFFFFF"),  // White
            TColor::GetColor("#FFF7BC"),  // Pale Yellow
            TColor::GetColor("#FEE391"),  // Light Yellow
            TColor::GetColor("#FEC44F"),  // Yellow
            TColor::GetColor("#FB9A29"),  // Orange
            TColor::GetColor("#EC7014"),  // Dark Orange
            TColor::GetColor("#CC4C02"),  // Brown
            TColor::GetColor("#993404"),  // Dark Brown
            TColor::GetColor("#662506")   // Very Dark Brown
        };
        colors = whOrBrColors;
    } else if (id == 10) {  // Iridescent
        cout << "Iridescent Palette!" << endl;
        nColors = 25;
        Int_t iridescentColors[] = {
            TColor::GetColor("#FEFBE9"),  // Very Pale Cream
            TColor::GetColor("#FCF7D5"),  // Pale Cream
            TColor::GetColor("#F5F3C1"),  // Light Cream
            TColor::GetColor("#EAF0B5"),  // Light Yellow
            TColor::GetColor("#DDECBF"),  // Pale Yellow
            TColor::GetColor("#D0E7CA"),  // Light Green
            TColor::GetColor("#B8E2C8"),  // Pale Green
            TColor::GetColor("#A3D6D1"),  // Light Teal
            TColor::GetColor("#8CC9D3"),  // Light Blue
            TColor::GetColor("#7FC2D7"),  // Blue
            TColor::GetColor("#6DACD7"),  // Medium Blue
            TColor::GetColor("#4E9ACD"),  // Dark Blue
            TColor::GetColor("#318CC0"),  // Darker Blue
            TColor::GetColor("#2F7C8C"),  // Teal
            TColor::GetColor("#3B4A6D"),  // Dark Teal
            TColor::GetColor("#504C77"),  // Grayish Blue
            TColor::GetColor("#7E7B8F"),  // Dark Gray
            TColor::GetColor("#9E8B98"),  // Grayish Pink
            TColor::GetColor("#BC7D7E"),  // Salmon
            TColor::GetColor("#E5A36D"),  // Light Orange
            TColor::GetColor("#F2B671"),  // Orange
            TColor::GetColor("#F1C2A8"),  // Peach
            TColor::GetColor("#F2A02A"),  // Yellow-Orange
            TColor::GetColor("#F08D22"),  // Dark Orange
            TColor::GetColor("#F06C1D"),  // Darker Orange
            TColor::GetColor("#CF3F2C")   // Red
        };
        colors = iridescentColors;
    } else {
        cout << "Unknown palette ID!" << endl;
        return;
    }
    
    gStyle->SetPalette(nColors, colors);
    gPad->Modified();
    gPad->Update();
}



void definePalette(){
    cout<<"\n;\n;\n--> To apply Palette: \n void SetWhichPalette(int id = 1 2 or 3).\n--> To test:\ndefinePalette_example() . "<<endl;
}

// Define the main function for ROOT to execute
void definePalette_example() {
    // Apply Wong's color palette
    SetWhichPalette(3);
    
    // Create a 2D histogram
    TH2F *h2 = new TH2F("h2", "Example 2D Histogram", 100, -4, 4, 100, -4, 4);

    // Fill 2D histogram with random data
    for (int i = 0; i < 10000; ++i) {
        h2->Fill(gRandom->Gaus(0, 1), gRandom->Gaus(0, 1));
    }

    // Create a canvas and draw the 2D histogram
    TCanvas *c = new TCanvas("c", "Canvas", 800, 600);
    h2->Draw("COLZ");  // Use COLZ option to show the color palette

    // Update the canvas to ensure it is drawn properly
    c->Update();
}
