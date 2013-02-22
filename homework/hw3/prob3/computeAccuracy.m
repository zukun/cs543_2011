function computeAccuracy(test_labels, predicted_labels)

disp(['Number of test examples: ' num2str(numel(test_labels))]);
disp(['Possible labels: ' num2str(unique(test_labels(:))')])
disp(['Accuracy = ' num2str(mean(test_labels(:)==predicted_labels(:)))]);
