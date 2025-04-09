#include <iostream>
#include <TDatime.h>
#include <TGraph.h>

Double_t calculateCharge(TGraph *gr, int month1, int day1, int hour1, int min1,
                     int month2, int day2, int hour2, int min2) {
    // Define the start and end datetimes

    int sec1 = 0;
    int sec2 = 0;
    TDatime datetime1(2023, month1, day1, hour1, min1, sec1);  // Start datetime (year fixed to 2023)
    TDatime datetime2(2023, month2, day2, hour2, min2, sec2);  // End datetime (year fixed to 2023)

    // Convert datetimes to timestamps
    Double_t timestamp1 = datetime1.Convert();
    Double_t timestamp2 = datetime2.Convert();

    // Initialize variables for total charge
    Double_t total_charge = 0.0;

    // Iterate through graph points to calculate charge
    for (Int_t i = 0; i < gr->GetN(); ++i) {
        Double_t x, y;
        gr->GetPoint(i, x, y);

        // Check if the point is within the desired interval
        if (x >= timestamp1 && x <= timestamp2) {
            // Calculate partial charge (current in µA multiplied by time in seconds)
            if(y>0)total_charge += y;
        }
    }

    // Display the total charge
    //std::cout << "Charge between " << datetime1.AsString() << " and " << datetime2.AsString() << ": " << total_charge << " µC" << std::endl;

    // Format and display the start and end datetimes
    auto formatDatetime = [](const TDatime& dt) {
        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(2) << dt.GetDay() << "/"
            << std::setfill('0') << std::setw(2) << dt.GetMonth() << " "
            << std::setfill('0') << std::setw(2) << dt.GetHour() << ":"
            << std::setfill('0') << std::setw(2) << dt.GetMinute();
        return oss.str();
    };

    std::cout << "Charge between " << formatDatetime(datetime1) << " and " << formatDatetime(datetime2)
              << ": " << total_charge << " µC" << std::endl;


    return total_charge;
}
