void clone_branch() {
    // Open the input file
    TFile* inputFile = TFile::Open("r0010_000a.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    // Get the tree from the file
    TTree* inputTree = dynamic_cast<TTree*>(inputFile->Get("AD"));
    if (!inputTree) {
        std::cerr << "TTree 'AD' not found in file!" << std::endl;
        inputFile->Close();
        return;
    }

    // Set up a new file to save the cloned branch
    TFile* outputFile = TFile::Open("cloned_branch.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output file!" << std::endl;
        inputFile->Close();
        return;
    }

    // Clone only the desired branch
    TTree* clonedTree = inputTree->CloneTree(0); // Create empty clone (structure only)
    inputTree->SetBranchStatus("*", 0); // Disable all branches
    inputTree->SetBranchStatus("Medley_1_dE2", 1); // Enable only the desired branch

    // Create a new tree with the selected branch
    clonedTree = inputTree->CloneTree(); // Clone selected entries

    // Write the tree to the output file
    outputFile->cd();
    clonedTree->Write("AD");

    // Clean up
    outputFile->Close();
    inputFile->Close();

    std::cout << "Branch 'Medley_1_dE2' cloned successfully to cloned_branch.root" << std::endl;
}